'''
This part of the workflow starts from files

  - pre-processed/sequences.fasta
  - pre-processed/metadata.tsv
  - pre-processed/sequence_index.tsv

and produces files

  - builds/{build_name}/sequences.fasta
  - builds/{build_name}/metadata.tsv

'''

localrules: freeze_archive_for_build, pango_update

build_dir = config.get("build_dir", "builds")

rule prepare_build:
    input:
        sequences = build_dir + "/{build_name}/sequences.fasta",
        metadata = build_dir + "/{build_name}/metadata.tsv"

rule freeze_archive_for_build:
    output:
        sequences = "freezed/pre-processed/filtered.fasta.xz",
        metadata = "freezed/pre-processed/metadata.tsv",
        sequence_index = "freezed/pre-processed/sequence_index.tsv",
    shell:
        """
            rm -rf freezed
            cp -r ../ncov-simple/archive freezed
        """

rule subsample:
    message:
        """
        Subsample all sequences by '{wildcards.subsample}' scheme for build '{wildcards.build_name}' with the following parameters:
        """
    input:
        sequences = "freezed/pre-processed/filtered.fasta.xz",
        metadata = "freezed/pre-processed/metadata.tsv",
        sequence_index = "freezed/pre-processed/sequence_index.tsv",
        include = config["files"]["include"],
    output:
        sequences = build_dir + "/{build_name}/sample-{subsample}.fasta",
        strains=build_dir + "/{build_name}/sample-{subsample}.txt",
    log:
        "logs/subsample_{build_name}_{subsample}.txt"
    benchmark:
        "benchmarks/subsample_{build_name}_{subsample}.txt"
    params:
        filter_arguments = lambda w: config["builds"][w.build_name]["subsamples"][w.subsample]['filters'],
        date = (datetime.date.today() + datetime.timedelta(days=1)).strftime("%Y-%m-%d"),
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024)

    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --sequence-index {input.sequence_index} \
            --include {input.include} \
            --max-date {params.date} \
            {params.filter_arguments} \
            --output {output.sequences} \
            --output-strains {output.strains} 2>&1 | tee {log}
        """

rule sequence_select:
    input:
        sequences = "freezed/pre-processed/filtered.fasta.xz",
        strains = lambda w: [build_dir + f"/{w.build_name}/sample-{subsample}.txt"
                   for subsample in config["builds"][w.build_name]["subsamples"]]
    output:
        sequences=build_dir + "/{build_name}/picked_pango.fasta",
        strains = build_dir + "/{build_name}/strains.txt"
    shell:
        """
        cat {input.strains} \
        | sort -u > {output.strains}
        xzcat {input.sequences} \
        | seqkit grep -f {output.strains} -o {output.sequences}
        """

rule combine_subsamples:
    # Similar to rule combine_input_metadata, this rule should only be run if multiple inputs are being used (i.e. multiple origins)
    message:
        """
        Combine and deduplicate aligned & filtered FASTAs from multiple origins in preparation for subsampling: {input}.
        """
    input: rules.sequence_select.output.sequences
    output:
        sequences = build_dir + "/{build_name}/sequences_raw.fasta",
    benchmark:
        "benchmarks/combine_subsamples_{build_name}.txt"
    shell:
        """
        python3 scripts/combine-and-dedup-fastas.py --input {input} --output {output}
        """

rule pango_assignments_default:
    input:
        sequences = rules.combine_subsamples.output.sequences,
    output:
        assignments = build_dir + "/{build_name}/pango_default.csv",
    log:
        "logs/pango_default_{build_name}.txt"
    shell:
        """
        pangolin --outfile {output.assignments} --analysis-mode pangolearn {input.sequences} 2>&1 | \
        tee {log}
        """

rule pango_assignments_usher:
    input:
        sequences = rules.combine_subsamples.output.sequences,
    output:
        assignments = build_dir + "/{build_name}/pango_usher.csv",
    log:
        "logs/pango_usher_{build_name}.txt"
    shell:
        """
        pangolin --outfile {output.assignments} {input.sequences} 2>&1 | \
        tee {log}
        """

rule extract_metadata:
    input:
        strains = lambda w: [build_dir + f"/{w.build_name}/sample-{subsample}.txt"
                   for subsample in config["builds"][w.build_name]["subsamples"]],
        metadata = "freezed/pre-processed/metadata.tsv"
    output:
        metadata = build_dir + "/{build_name}/metadata_pre_pango.tsv",
    params:
        adjust = lambda w: config["builds"][w.build_name].get("metadata_adjustments",{}),
    benchmark:
        "benchmarks/extract_metadata_{build_name}.txt"
    run:
        import pandas as pd
        strains = set()
        for f in input.strains:
            with open(f) as fh:
                strains.update([x.strip() for x in fh if x[0]!='#'])

        d = pd.read_csv(input.metadata, index_col='strain', sep='\t').loc[list(strains)]
        if len(params.adjust):
            for adjustment  in params.adjust:
                ind = d.eval(adjustment["query"])
                d.loc[ind, adjustment['dst']] = d.loc[ind, adjustment['src']]

        d.to_csv(output.metadata, sep='\t')

rule join_pangolins:
    input:
        default = rules.pango_assignments_default.output.assignments,
        usher = rules.pango_assignments_usher.output.assignments,
        metadata = rules.extract_metadata.output.metadata
    output:
        metadata = rules.prepare_build.input.metadata
    run:
        import pandas as pd
        raw_meta = pd.read_csv(input.metadata, index_col='strain', sep='\t')
        default = pd.read_csv(input.default)
        default.set_index('taxon',inplace=True)
        default.rename(columns={'lineage':'pango_default'}, inplace=True)
        meta = raw_meta.join(default['pango_default'], how='left')
        usher = pd.read_csv(input.usher)
        usher.set_index('taxon',inplace=True)
        usher.rename(columns={'lineage':'pango_usher'}, inplace=True)
        meta = meta.join(usher['pango_usher'], how='left')
        meta.to_csv(output.metadata, sep='\t')

rule exclude_outliers:
    input:
        sequences = rules.combine_subsamples.output.sequences,
        metadata = rules.prepare_build.input.metadata,
        exclude = "profiles/exclude.txt",
    output:
        sampled_sequences = rules.prepare_build.input.sequences,
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output {output.sampled_sequences}
        """
