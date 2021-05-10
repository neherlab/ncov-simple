'''
This part of the workflow downloads files from S3

  - "data/{origin}/sequences.fasta.gz"
  - "data/{origin}/metadata.tsv"

and produces

  - pre-processed/sequences.fasta
  - pre-processed/metadata.tsv

'''

import os
from snakemake.remote.S3 import RemoteProvider as S3Provider
S3 = S3Provider()

localrules: download_sequences, download_metadata, download_exclude

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
    message: "Downloading sequences from {input[0]} -> {output[0]}"
    input:
        lambda w: S3.remote(config['origins'][w.origin]['sequences'])
    output:
        "data/{origin}/sequences.fasta.gz"
    shell: "mv {input} {output}"

rule download_metadata:
    message: "Downloading metadata from {input} -> {output}"
    input:
        lambda w: S3.remote(config['origins'][w.origin]['metadata'])
    output:
        metadata = "data/{origin}/metadata.tsv"
    params:
        deflate = lambda w: _infer_decompression(config['origins'][w.origin]['metadata'])
    shell: "{params.deflate} {input} > {output:q}"

rule download_exclude:
    message: "Downloading exclude from {input} -> {output}"
    output:
        "data/{origin}/exclude.txt"
    params:
        source = lambda w: config["origins"][w.origin]['exclude']
    shell: "curl {params.source} -o {output}"

rule prealign:
    message:
        """
        Aligning sequences to {input.reference}
            - gaps relative to reference are considered real
        """
    input:
        sequences = "data/{origin}/sequences.fasta.gz",
        genemap = config["files"]["annotation"],
        reference = config["files"]["alignment_reference"]
    output:
        alignment = "pre-processed/{origin}/alignment.fasta.xz",
        insertions = "pre-processed/{origin}/insertions.tsv",
        translations = expand("pre-processed/{{origin}}/translations/seqs.gene.{gene}.fasta.xz", gene=config.get('genes', ['S']))
    params:
        outdir = "pre-processed/{origin}/translations",
        genes = ','.join(config.get('genes', ['S'])),
        basename = "seqs",
        deflate = lambda w: _infer_decompression(".gz")
    log:
        "logs/prealign_{origin}.txt"
    benchmark:
        "benchmarks/align_{origin}.txt"
    conda: config["conda_environment"]
    threads: 8
    resources:
        mem_mb=3000
    shell:
        """
        {params.deflate} {input.sequences} | nextalign \
            --jobs={threads} \
            --reference {input.reference} \
            --genemap {input.genemap} \
            --genes {params.genes} \
            --sequences /dev/stdin \
            --output-dir {params.outdir} \
            --output-basename {params.basename} \
            --output-fasta /dev/stdout \
            --output-insertions {output.insertions} | xz -2c - > {output.alignment} &&\
            xz -2 {params.outdir}/*fasta
        """

rule mask:
    message:
        """
        Mask bases in alignment {input.alignment}
          - masking {params.mask_arguments}
        """
    input:
        alignment = "pre-processed/{origin}/alignment.fasta.xz"
    output:
        alignment = "pre-processed/{origin}/masked.fasta.xz"
    log:
        "logs/mask_{origin}.txt"
    benchmark:
        "benchmarks/mask_{origin}.txt"
    params:
        mask_arguments = lambda w: config["origins"][w.origin].get("mask","")
    conda: config["conda_environment"]
    shell:
        """
        augur mask \
            --sequences {input.alignment} \
            {params.mask_arguments} \
            --output {output.alignment} 2>&1 | tee {log}
        """

rule filter:
    message:
        """
        Filtering alignment {input.sequences} -> {output.sequences}
          - excluding strains in {input.exclude}
          - including strains in {input.include}
        """
    input:
        sequences = "pre-processed/{origin}/masked.fasta.xz",
        metadata = "data/{origin}/metadata.tsv",
        include = "defaults/include.txt",
        exclude = "data/{origin}/exclude.txt"
    output:
        sequences = "pre-processed/{origin}/filtered.fasta.xz"
    log:
        "logs/filtered{origin}.txt"
    benchmark:
        "benchmarks/filter{origin}.txt"
    params:
        filter_arguments = lambda w: config["origins"][w.origin].get("filters","")
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include {input.include} \
            --exclude {input.exclude} \
            {params.filter_arguments} \
            --output {output.sequences} 2>&1 | tee {log}
        """


rule combine_bulk_sequences:
    input:
        [f"pre-processed/{origin}/filtered.fasta.xz" for origin in config["origins"]]
    output:
        "pre-processed/filtered.fasta.xz"
    run:
        if len(input)==1:
            shell(f"cp {input} {output}")

rule combine_bulk_metadata:
    input:
        [f"data/{origin}/metadata.tsv" for origin in config["origins"]]
    output:
        "pre-processed/metadata.tsv"
    run:
        if len(input)==1:
            shell(f"cp {input} {output}")

rule index_sequences:
    message:
        """
        Index sequence composition for faster filtering.
        """
    input:
        sequences = rules.combine_bulk_sequences.output
    output:
        sequence_index = "pre-processed/sequence_index.tsv"
    log:
        "logs/index_sequences.txt"
    benchmark:
        "benchmarks/index_sequences.txt"
    conda: config["conda_environment"]
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index} 2>&1 | tee {log}
        """
