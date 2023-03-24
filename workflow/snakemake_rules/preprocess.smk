"""
This part of the workflow downloads files from S3

  - "data/{origin}/sequences.fasta.xz"
  - "data/{origin}/metadata.tsv"

and produces

  - pre-processed/filtered.fasta.xz
  - pre-processed/metadata.tsv

"""

import os


localrules:
    download_clade_emergence_dates,
    download_sequences,
    download_mutation_summary,
    download_metadata,
    download_exclude,
    download_clades,
    preprocess,
    download_lat_longs,
    download_color_ordering,
    download_mutational_fitness_map,


rule preprocess:
    input:
        sequences="pre-processed/filtered.fasta.xz",
        metadata="pre-processed/metadata.tsv",
        sequence_index="pre-processed/sequence_index.tsv",
    params:
        slack_hook=config.get("slackHook", "google.com"),
    shell:
        """
        rm -rf archive/pre-processed
        cp -r pre-processed archive/pre-processed
        curl -X POST -H 'Content-type: application/json' \
        --data '{{"text":"Preprocessing done"}}' \
        {params.slack_hook}
        """


def _infer_decompression(input):
    """
    Returns a shell command to decompress the piped stream,
    which will itself produce a stream of decompressed data to stdout.
    If no decompression is needed, returns `cat`.
    NOTE: a lot of this will become unnecessary once `augur` handles
    compressed sequence inputs.
    """
    if input.endswith(".xz"):
        return "xz -dcq"
    if input.endswith(".gz"):
        return "gunzip -cq"
    return "cat"


rule download_sequences:
    message:
        "Downloading sequences from {params.address} -> {output[0]}"
    params:
        address=lambda w: config["origins"][w.origin]["sequences"],
    output:
        "pre-processed/{origin}/alignment.fasta.xz",
    shell:
        "aws s3 cp {params.address} {output}"


rule download_metadata:
    message:
        "Downloading metadata from {params.address} -> {output}"
    params:
        deflate=lambda w: _infer_decompression(config["origins"][w.origin]["metadata"]),
        address=lambda w: config["origins"][w.origin]["metadata"],
    output:
        metadata="data/{origin}/metadata.tsv",
    shell:
        "aws s3 cp {params.address} - | {params.deflate} {input} > {output:q}"


rule download_mutation_summary:
    message:
        "Downloading mutation summary from {params.address} -> {output}"
    params:
        deflate=lambda w: _infer_decompression(
            config["origins"][w.origin]["mutation_summary"]
        ),
        address=lambda w: config["origins"][w.origin]["mutation_summary"],
    output:
        mutation_summary="pre-processed/{origin}/mutation_summary.tsv",
    shell:
        "aws s3 cp {params.address} - | {params.deflate} {input} > {output:q}"


rule download_exclude:
    message:
        "Downloading exclude from {params.source} -> {output}"
    output:
        "data/{origin}/exclude.txt",
    params:
        source=lambda w: config["origins"][w.origin]["exclude"],
    shell:
        "curl {params.source} -o {output}"


rule download_clades:
    message:
        "Downloading clade definitions from {params.source} -> {output}"
    output:
        config["files"]["clades"],
    params:
        source=config["data_source"]["clades"],
    shell:
        "curl {params.source} -o {output}"


rule download_color_ordering:
    message:
        "Downloading clade definitions from {params.source} -> {output}"
    output:
        config["files"]["ordering"],
    params:
        source=config["data_source"]["color_ordering"],
    shell:
        "curl {params.source} -o {output}"


rule download_lat_longs:
    message:
        "Downloading clade definitions from {params.source} -> {output}"
    output:
        config["files"]["lat_longs"],
    params:
        source=config["data_source"]["lat_longs"],
    shell:
        "curl {params.source} -o {output}"


rule download_mutational_fitness_map:
    message:
        "Downloading mutational fitness map from {params.source} -> {output}"
    output:
        config["files"]["mut_fit"],
    params:
        source=config["data_source"]["mut_fit"],
    shell:
        "curl {params.source} -o {output}"


rule download_clade_emergence_dates:
    message:
        "Downloading clade emergence dates from {params.source} -> {output}"
    output:
        config["files"]["clade_emergence_dates"],
    params:
        source=config["data_source"]["clade_emergence_dates"],
    shell:
        "curl {params.source} -o {output}"


rule diagnostic:
    message:
        "Scanning metadata {input.metadata} for problematic sequences. Removing sequences with >{params.clock_filter} deviation from the clock and with more than {params.snp_clusters}."
    input:
        metadata="data/{origin}/metadata.tsv",
        clade_emergence_dates=rules.download_clade_emergence_dates.output,
    output:
        to_exclude="pre-processed/{origin}/problematic_exclude.txt",
        exclude_reasons="pre-processed/{origin}/exclude_reasons.txt",
    params:
        clock_filter=12,
        snp_clusters=1,
        contamination=5,
        clock_plus_rare=45,
    log:
        "logs/diagnostics_{origin}.txt",
    benchmark:
        "benchmarks/diagnostics_{origin}.txt"
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=12000,
    shell:
        """
        python3 scripts/diagnostic.py \
            --metadata {input.metadata} \
            --clock-filter {params.clock_filter} \
            --contamination {params.contamination} \
            --clade-emergence-dates {input.clade_emergence_dates} \
            --snp-clusters {params.snp_clusters} \
            --output-exclusion-list {output.to_exclude} \
            --output-exclusion-reasons {output.exclude_reasons} \
            2>&1 | tee {log}
        """


rule filter:
    input:
        metadata="data/{origin}/metadata.tsv",
        include="defaults/include.txt",
        exclude="data/{origin}/exclude.txt",
        problematic="pre-processed/{origin}/problematic_exclude.txt",
    output:
        metadata="pre-processed/{origin}/metadata.tsv",
    log:
        "logs/filtered{origin}.txt",
    benchmark:
        "benchmarks/filter{origin}.txt"
    params:
        filter_arguments=lambda w: config["origins"][w.origin].get("filters", ""),
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024),
    shell:
        """
        augur filter \
            --metadata {input.metadata} \
            --include {input.include} \
            --exclude {input.exclude} {input.problematic} \
            {params.filter_arguments} \
            --output-metadata {output.metadata} 2>&1 | tee {log};
        """


rule combine_bulk_sequences:
    input:
        sequences=[
            f"pre-processed/{origin}/alignment.fasta.xz"
            for origin in config["origins"]
        ],
        # mutation_summary = [f"pre-processed/{origin}/mutation_summary.tsv" for origin in config["origins"]]
    output:
        rules.preprocess.input.sequences,
    shell:
        """
        cp {input.sequences} {output}
        """


rule combine_bulk_metadata:
    input:
        [f"pre-processed/{origin}/metadata.tsv" for origin in config["origins"]],
    output:
        rules.preprocess.input.metadata,
    shell:
        """
        cp {input} {output}
        """


rule split_sequences:
    input:
        sequences="pre-processed/filtered.fasta.xz",
    output:
        sequences=expand(
            "pre-processed/split/filtered.part_00{part_no}.fasta.zst",
            part_no=range(1, 10),
        ),
    log:
        "logs/split_sequences.txt",
    shell:
        """
        seqkit split2 \
            {input.sequences} \
            -p 9 \
            -f \
            -e .zst \
            -O pre-processed/split
        """


rule index_sequences:
    """
    Index sequence composition for faster filtering.
    """
    input:
        sequences="pre-processed/split/filtered.part_00{part_no}.fasta.zst",
    output:
        sequence_index="pre-processed/index_00{part_no}.tsv",
    log:
        "logs/index_sequences_{part_no}.txt",
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index} 2>&1 | tee {log}
        """


rule combine_index:
    """
    Combine sequence composition indices.
    """
    input:
        sequence_index=expand(
            rules.index_sequences.output.sequence_index, part_no=range(1, 10)
        ),
    output:
        sequence_index="pre-processed/sequence_index.tsv",
    shell:
        """
        keep-header {input.sequence_index} -- cat > {output.sequence_index}
        """
