'''
This part of the workflow starts from files

  - pre-processed/sequences.fasta
  - pre-processed/metadata.tsv

and produces files

  - builds/{build_name}/sequences.fasta
  - builds/{build_name}/metadata.tsv

'''
rule prepare_build:
    input:
        sequences = "builds/{build_name}/sequences.fasta",
        metadata = "builds/{build_name}/metadata.tsv"

def _get_priority_file(w):
    if "priorities" in config["builds"][w.build_name]["subsamples"][w.subsample]:
        return f"builds/{w.build_name}/priorities_{config['builds'][w.build_name]['subsamples'][w.subsample].get('priorities')}.tsv"
    else:
        return []

def _get_priority_argument(w):
    f = _get_priority_file(w)
    if f:
        return "--priority " + f
    else:
        return ""

rule subsample:
    message:
        """
        Subsample all sequences by '{wildcards.subsample}' scheme for build '{wildcards.build_name}' with the following parameters:
        """
    input:
        sequences = "pre-processed/filtered.fasta.xz",
        metadata = "pre-processed/metadata.tsv",
        sequence_index = "pre-processed/sequence_index.tsv",
        include = config["files"]["include"],
        priorities = _get_priority_file
    output:
        sequences = "builds/{build_name}/sample-{subsample}.fasta",
        strains="builds/{build_name}/sample-{subsample}.txt",
    log:
        "logs/subsample_{build_name}_{subsample}.txt"
    benchmark:
        "benchmarks/subsample_{build_name}_{subsample}.txt"
    params:
        filter_arguments = lambda w: config["builds"][w.build_name]["subsamples"][w.subsample]['filters'],
        priorities = _get_priority_argument
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --sequence-index {input.sequence_index} \
            --include {input.include} \
            {params.filter_arguments} \
            {params.priorities} \
            --output {output.sequences} \
            --output-strains {output.strains} 2>&1 | tee {log}
        """

rule proximity_score:
    message:
        """
        determine priority for inclusion in as phylogenetic context by
        genetic similiarity to sequences in focal set for build '{wildcards.build_name}'.
        """
    input:
        alignment = "pre-processed/filtered.fasta.xz",
        reference = config["files"]["alignment_reference"],
        focal_alignment = "builds/{build_name}/sample-{focus}.fasta"
    output:
        proximities = "builds/{build_name}/proximity_{focus}.tsv"
    log:
        "logs/subsampling_proximity_{build_name}_{focus}.txt"
    benchmark:
        "benchmarks/proximity_score_{build_name}_{focus}.txt"
    params:
        chunk_size=10000,
        ignore_seqs = config['refine']['root']
    resources:
        # Memory scales at ~0.15 MB * chunk_size (e.g., 0.15 MB * 10000 = 1.5GB).
        mem_mb=4000
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/get_distance_to_focal_set.py \
            --reference {input.reference} \
            --alignment {input.alignment} \
            --focal-alignment {input.focal_alignment} \
            --ignore-seqs {params.ignore_seqs} \
            --chunk-size {params.chunk_size} \
            --output {output.proximities} 2>&1 | tee {log}
        """

rule priority_score:
    input:
        proximity = rules.proximity_score.output.proximities,
        sequence_index = rules.index_sequences.output.sequence_index,
    output:
        priorities = "builds/{build_name}/priorities_{focus}.tsv"
    benchmark:
        "benchmarks/priority_score_{build_name}_{focus}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/priorities.py \
            --sequence-index {input.sequence_index} \
            --proximities {input.proximity} \
            --output {output.priorities} 2>&1 | tee {log}
        """


rule combine_subsamples:
    # Similar to rule combine_input_metadata, this rule should only be run if multiple inputs are being used (i.e. multiple origins)
    message:
        """
        Combine and deduplicate aligned & filtered FASTAs from multiple origins in preparation for subsampling: {input}.
        """
    input:
        lambda w: [f"builds/{w.build_name}/sample-{subsample}.fasta"
                   for subsample in config["builds"][w.build_name]["subsamples"]]
    output:
        rules.prepare_build.input.sequences
    benchmark:
        "benchmarks/combine_subsamples_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/combine-and-dedup-fastas.py --input {input} --output {output}
        """


rule extract_metadata:
    input:
        strains = lambda w: [f"builds/{w.build_name}/sample-{subsample}.txt"
                   for subsample in config["builds"][w.build_name]["subsamples"]],
        metadata = "pre-processed/metadata.tsv"
    output:
        rules.prepare_build.input.metadata
    benchmark:
        "benchmarks/extract_metadata_{build_name}.txt"
    run:
        import pandas as pd
        strains = set()
        for f in input.strains:
            with open(f) as fh:
                strains.update([x.strip() for x in fh if x[0]!='#'])

        pd.read_csv(input.metadata, index_col=0, sep='\t').loc[list(strains)].to_csv(output[0], sep='\t')


