"""
This part of the workflow starts from files

  - builds/{build_name}/sequences.fasta
  - builds/{build_name}/metadata.tsv

and produces files

  - auspice/ncov_{build_name}.json
  - auspice/ncov_{build_name}-tip-frequencies.json
  - auspice/ncov_{build_name}-root-sequence.json

"""


localrules:
    timestamped_build,
    include_hcov19_prefix,


build_dir = config.get("build_dir", "builds")
auspice_dir = config.get("auspice_dir", "auspice")
auspice_prefix = config.get("auspice_prefix", "ncov")


rule align:
    input:
        sequences=build_dir + "/{build_name}/sequences.fasta",
        genemap=config["files"]["annotation"],
        reference=config["files"]["alignment_reference"],
    output:
        alignment=build_dir + "/{build_name}/aligned.fasta",
        insertions=build_dir + "/{build_name}/insertions.tsv",
        translations=expand(
            build_dir + "/{{build_name}}/translations/aligned.gene.{gene}.fasta",
            gene=config.get("genes", ["S"]),
        ),
    params:
        outdir=lambda w: build_dir
        + f"/{w.build_name}/"
        + "translations/aligned.gene.{gene}.fasta",
        genes=",".join(config.get("genes", ["S"])),
    log:
        "logs/align_{build_name}.txt",
    benchmark:
        "benchmarks/align_{build_name}.txt"
    threads: 4
    resources:
        mem_mb=3000,
    shell:
        """
        nextalign run \
            --jobs={threads} \
            --input-ref {input.reference} \
            --input-gene-map {input.genemap} \
            --genes {params.genes} \
            {input.sequences} \
            --output-translations {params.outdir} \
            --output-fasta {output.alignment} \
            --output-insertions {output.insertions}
            > {log} 2>&1
        """


rule mask:
    input:
        alignment=rules.align.output.alignment,
    output:
        alignment=build_dir + "/{build_name}/masked.fasta",
    log:
        "logs/mask_{build_name}.txt",
    benchmark:
        "benchmarks/mask_{build_name}.txt"
    params:
        mask_arguments=lambda w: config.get("mask", ""),
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            {params.mask_arguments} \
            --output {output.alignment} 2>&1 | tee {log}
        """


rule tree:
    input:
        alignment=rules.mask.output.alignment,
    output:
        tree=build_dir + "/{build_name}/tree_raw.nwk",
    params:
        args=lambda w: config["tree"].get("tree-builder-args", "")
        if "tree" in config
        else "",
        exclude_sites=lambda w: f"--exclude-sites {config['files']['sites_to_mask']}"
        if "sites_to_mask" in config["files"]
        else "",
    log:
        "logs/tree_{build_name}.txt",
    benchmark:
        "benchmarks/tree_{build_name}.txt"
    threads: 8
    resources:
        # Multiple sequence alignments can use up to 40 times their disk size in
        # memory, especially for larger alignments.
        # Note that Snakemake >5.10.0 supports input.size_mb to avoid converting from bytes to MB.
        mem_mb=lambda wildcards, input: 40 * int(input.size / 1024 / 1024),
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --tree-builder-args {params.args} \
            --output {output.tree} \
            {params.exclude_sites} \
            --nthreads {threads} 2>&1 | tee {log}
        """


rule refine:
    input:
        tree=rules.tree.output.tree,
        alignment=rules.mask.output.alignment,
        metadata=build_dir + "/{build_name}/metadata.tsv",
    output:
        tree=build_dir + "/{build_name}/tree.nwk",
        node_data=build_dir + "/{build_name}/branch_lengths.json",
    log:
        "logs/refine_{build_name}.txt",
    benchmark:
        "benchmarks/refine_{build_name}.txt"
    threads: 1
    resources:
        # Multiple sequence alignments can use up to 15 times their disk size in
        # memory.
        # Note that Snakemake >5.10.0 supports input.size_mb to avoid converting from bytes to MB.
        mem_mb=lambda wildcards, input: 15 * int(input.size / 1024 / 1024),
    params:
        root=config["refine"]["root"],
        clock_rate=config["refine"].get("clock_rate", 0.0007),
        clock_std_dev=config["refine"].get("clock_std_dev", 0.003),
        coalescent=config["refine"].get("coalescent", "opt"),
        date_inference=config["refine"].get("date_inference", "marginal"),
        divergence_unit=config["refine"].get("divergence_unit", "mutations"),
        clock_filter_iqd=config["refine"].get("clock_filter_iqd", 4),
        keep_polytomies="--keep-polytomies"
        if config["refine"].get("keep_polytomies", False)
        else "",
        timetree="" if config["refine"].get("no_timetree", False) else "--timetree",
    shell:
        """
        augur refine \
            --use-fft \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --root {params.root} \
            {params.timetree} \
            {params.keep_polytomies} \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --coalescent {params.coalescent} \
            --date-inference {params.date_inference} \
            --divergence-unit {params.divergence_unit} \
            --date-confidence \
            --no-covariance \
            --clock-filter-iqd {params.clock_filter_iqd} 2>&1 | tee {log}
        """


rule ancestral:
    input:
        tree=rules.refine.output.tree,
        alignment=rules.align.output.alignment,
    output:
        node_data=build_dir + "/{build_name}/nt_muts.json",
    log:
        "logs/ancestral_{build_name}.txt",
    benchmark:
        "benchmarks/ancestral_{build_name}.txt"
    params:
        inference=config["ancestral"]["inference"],
    resources:
        # Multiple sequence alignments can use up to 15 times their disk size in
        # memory.
        # Note that Snakemake >5.10.0 supports input.size_mb to avoid converting from bytes to MB.
        mem_mb=lambda wildcards, input: 15 * int(input.size / 1024 / 1024),
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} \
            --infer-ambiguous 2>&1 | tee {log}
        """


rule translate:
    input:
        tree=rules.refine.output.tree,
        translations=lambda w: rules.align.output.translations,
        reference=config["files"]["alignment_reference"],
        genemap=config["files"]["annotation"],
    output:
        node_data=build_dir + "/{build_name}/aa_muts.json",
        translations=expand(
            build_dir
            + "/{{build_name}}/translations/aligned.gene.{gene}_withInternalNodes.fasta",
            gene=config.get("genes", ["S"]),
        ),
    params:
        genes=config.get("genes", "S"),
    log:
        "logs/aamuts_{build_name}.txt",
    benchmark:
        "benchmarks/aamuts_{build_name}.txt"
    resources:
        # Multiple sequence alignments can use up to 15 times their disk size in
        # memory.
        # Note that Snakemake >5.10.0 supports input.size_mb to avoid converting from bytes to MB.
        mem_mb=lambda wildcards, input: 15 * int(input.size / 1024 / 1024),
    shell:
        """
        python3 scripts/explicit_translation.py \
            --tree {input.tree} \
            --annotation {input.genemap} \
            --reference {input.reference} \
            --translations {input.translations:q} \
            --genes {params.genes} \
            --output {output.node_data} 2>&1 | tee {log}
        """


rule traits:
    input:
        tree=rules.refine.output.tree,
        metadata=build_dir + "/{build_name}/metadata.tsv",
    output:
        node_data=build_dir + "/{build_name}/traits.json",
    log:
        "logs/traits_{build_name}.txt",
    benchmark:
        "benchmarks/traits_{build_name}.txt"
    params:
        columns=config["traits"]["columns"],
        sampling_bias_correction=config["traits"]["sampling_bias_correction"],
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024),
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence \
            --sampling-bias-correction {params.sampling_bias_correction} 2>&1 | tee {log}
        """


rule clades:
    input:
        tree=rules.refine.output.tree,
        aa_muts=rules.translate.output.node_data,
        nuc_muts=rules.ancestral.output.node_data,
        clades=config["files"]["clades"],
    output:
        node_data=build_dir + "/{build_name}/clades.json",
    log:
        "logs/clades_{build_name}.txt",
    benchmark:
        "benchmarks/clades_{build_name}.txt"
    resources:
        # Memory use scales primarily with size of the node data.
        mem_mb=lambda wildcards, input: 3 * int(input.size / 1024 / 1024),
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """


rule tip_frequencies:
    input:
        tree=rules.refine.output.tree,
        metadata=build_dir + "/{build_name}/metadata.tsv",
    output:
        tip_frequencies_json=build_dir + "/{build_name}/tip-frequencies.json",
    log:
        "logs/tip_frequencies_{build_name}.txt",
    benchmark:
        "benchmarks/tip_frequencies_{build_name}.txt"
    params:
        min_date=config["frequencies"]["min_date"],
        max_date=lambda w: datetime.datetime.today().strftime("%Y-%m-%d"),
        pivot_interval=config["frequencies"]["pivot_interval"],
        pivot_interval_units=config["frequencies"]["pivot_interval_units"],
        narrow_bandwidth=config["frequencies"]["narrow_bandwidth"],
        proportion_wide=config["frequencies"]["proportion_wide"],
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024),
    shell:
        """
        augur frequencies \
            --method kde \
            --metadata {input.metadata} \
            --tree {input.tree} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --pivot-interval {params.pivot_interval} \
            --pivot-interval-units {params.pivot_interval_units} \
            --narrow-bandwidth {params.narrow_bandwidth} \
            --proportion-wide {params.proportion_wide} \
            --output {output.tip_frequencies_json} 2>&1 | tee {log}
        """


if "distances" in config:

    rule distances:
        input:
            tree=rules.refine.output.tree,
            alignment=build_dir
            + "/{build_name}/translations/aligned.gene.S_withInternalNodes.fasta",
            distance_maps=config["distances"]["maps"],
        params:
            genes="S",
            comparisons=config["distances"]["comparisons"],
            attribute_names=config["distances"]["attributes"],
        output:
            node_data=build_dir + "/{build_name}/distances.json",
        shell:
            """
            augur distance \
                --tree {input.tree} \
                --alignment {input.alignment} \
                --gene-names {params.genes} \
                --compare-to {params.comparisons} \
                --attribute-name {params.attribute_names} \
                --map {input.distance_maps} \
                --output {output}
            """


rule mutational_fitness:
    input:
        tree=rules.refine.output.tree,
        alignment=lambda w: rules.translate.output.translations,
        distance_maps=rules.download_mutational_fitness_map.output,
    output:
        node_data=build_dir + "/{build_name}/mutational_fitness.json",
    benchmark:
        "benchmarks/mutational_fitness_{build_name}.txt"
    log:
        "logs/mutational_fitness_{build_name}.txt",
    params:
        genes=" ".join(config.get("genes", ["S"])),
        compare_to="root",
        attribute_name="mutational_fitness",
    resources:
        mem_mb=2000,
    shell:
        """
        augur distance \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --gene-names {params.genes} \
            --compare-to {params.compare_to} \
            --attribute-name {params.attribute_name} \
            --map {input.distance_maps} \
            --output {output} 2>&1 | tee {log}
        """


rule colors:
    input:
        ordering=config["files"]["ordering"],
        color_schemes=config["files"]["color_schemes"],
        metadata=build_dir + "/{build_name}/metadata.tsv",
    output:
        colors=build_dir + "/{build_name}/colors.tsv",
    log:
        "logs/colors_{build_name}.txt",
    benchmark:
        "benchmarks/colors_{build_name}.txt"
    resources:
        # Memory use scales primarily with the size of the metadata file.
        # Compared to other rules, this rule loads metadata as a pandas
        # DataFrame instead of a dictionary, so it uses much less memory.
        mem_mb=lambda wildcards, input: 5 * int(input.metadata.size / 1024 / 1024),
    shell:
        """
        python3 scripts/assign-colors.py \
            --ordering {input.ordering} \
            --color-schemes {input.color_schemes} \
            --output {output.colors} \
            --metadata {input.metadata} 2>&1 | tee {log}
        """


rule recency:
    input:
        metadata=build_dir + "/{build_name}/metadata.tsv",
    output:
        node_data=build_dir + "/{build_name}/recency.json",
    log:
        "logs/recency_{build_name}.txt",
    benchmark:
        "benchmarks/recency_{build_name}.txt"
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=12000,
    shell:
        """
        python3 scripts/construct-recency-from-submission-date.py \
            --metadata {input.metadata} \
            --output {output} 2>&1 | tee {log}
        """


def _get_node_data_by_wildcards(wildcards):
    """Return a list of node data files to include for a given build's wildcards."""
    # Define inputs shared by all builds.
    wildcards_dict = dict(wildcards)
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.clades.output.node_data,
        rules.traits.output.node_data,
        rules.mutational_fitness.output.node_data,
        # rules.recency.output.node_data #Doesn't work anymore with augur.read_metadata deprecated
    ]
    if "distances" in config:
        inputs.append(rules.distances.output.node_data)

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards_dict) for input_file in inputs]
    return inputs


# Make switzerland build have swiss acknowledgments
if "CH-geneva" in config["builds"]:
    config["builds"]["switzerland"]["description"] = config["builds"]["CH-geneva"][
        "description"
    ]
    config["builds"]["switzerland"]["auspice_config"] = config["builds"]["CH-geneva"][
        "auspice_config"
    ]


rule export:
    input:
        tree=rules.refine.output.tree,
        metadata=build_dir + "/{build_name}/metadata.tsv",
        node_data=_get_node_data_by_wildcards,
        auspice_config=lambda w: config["builds"][w.build_name]["auspice_config"]
        if "auspice_config" in config["builds"][w.build_name]
        else config["files"]["auspice_config"],
        colors=lambda w: config["builds"][w.build_name]["colors"]
        if "colors" in config["builds"][w.build_name]
        else (
            config["files"]["colors"]
            if "colors" in config["files"]
            else rules.colors.output.colors.format(**w)
        ),
        lat_longs=config["files"]["lat_longs"],
        description=lambda w: config["builds"][w.build_name]["description"]
        if "description" in config["builds"][w.build_name]
        else config["files"]["description"],
        tip_freq_json=rules.tip_frequencies.output.tip_frequencies_json,
    output:
        auspice_json=auspice_dir + f"/{{build_name}}/raw_nohcov.json",
        root_sequence_json=auspice_dir
        + f"/{{build_name}}/raw_nohcov_root-sequence.json",
        tip_freq_json=auspice_dir + f"/{{build_name}}/raw_nohcov_tip-frequencies.json",
    log:
        "logs/export_{build_name}.txt",
    benchmark:
        "benchmarks/export_{build_name}.txt"
    params:
        title=lambda w: config["builds"][w.build_name].get(
            "title", "SARS-CoV-2 phylogeny"
        ),
    # resources:
    # Memory use scales primarily with the size of the metadata file.
    # mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024)

    wildcard_constraints:
        build_name="[^_]+(_[^_]+)?",
    shell:
        """
        export AUGUR_RECURSION_LIMIT=10000;
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --include-root-sequence \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --title {params.title:q} \
            --description {input.description} \
            --skip-validation \
            --output {output.auspice_json} 2>&1 | tee {log};
            cp {input.tip_freq_json} {output.tip_freq_json}
        """


rule infer_insertions:
    input:
        metadata=build_dir + "/{build_name}/metadata.tsv",
        tree=rules.refine.output.tree,
    output:
        build_dir + "/{build_name}/insertions.json",
    log:
        "logs/infer_insertions_{build_name}.txt",
    shell:
        """
        python3 scripts/reconstruct_insertions.py \
            --metadata {input.metadata} \
            --tree {input.tree} \
            --output {output} 2>&1 | tee {log}
        """


rule add_branch_labels:
    input:
        auspice_json=auspice_dir + f"/{{build_name}}/raw_nohcov.json",
        mutations=rules.translate.output.node_data,
        insertions=build_dir + "/{build_name}/insertions.json",
    output:
        auspice_json=auspice_dir + f"/{{build_name}}/nohcov_auspice.json",
    log:
        "logs/add_branch_labels_{build_name}.txt",
    wildcard_constraints:
        build_name="[^_]+(_[^_]+)?",
    shell:
        """
        python3 scripts/add_branch_labels.py \
            --input {input.auspice_json} \
            --mutations {input.mutations} \
            --insertions {input.insertions} \
            --output {output.auspice_json} 
        """


rule include_hcov19_prefix:
    input:
        auspice_json=auspice_dir + f"/{{build_name}}/nohcov_auspice.json",
        root_sequence_json=auspice_dir
        + f"/{{build_name}}/raw_nohcov_root-sequence.json",
        tip_freq_json=auspice_dir + f"/{{build_name}}/raw_nohcov_tip-frequencies.json",
    output:
        auspice_json=auspice_dir
        + f"/{{build_name}}/latest/{auspice_prefix}_{{build_name}}.json",
        tip_freq_json=auspice_dir
        + f"/{{build_name}}/latest/{auspice_prefix}_{{build_name}}_tip-frequencies.json",
        root_sequence_json=auspice_dir
        + f"/{{build_name}}/latest/{auspice_prefix}_{{build_name}}_root-sequence.json",
    log:
        "logs/include_hcov19_prefix_{build_name}.txt",
    wildcard_constraints:
        build_name="[^_]+(_[^_]+)?",
    shell:
        """
        python3 ./scripts/include_prefix.py \
            --input-auspice {input.auspice_json} \
            --input-tip-frequencies {input.tip_freq_json} \
            --prefix "hCoV-19/" \
            --output-auspice {output.auspice_json} \
            --output-tip-frequencies {output.tip_freq_json}; \
        cp {input.root_sequence_json} {output.root_sequence_json}; \
        rm {input.auspice_json} {input.tip_freq_json} {input.root_sequence_json}
    """


rule timestamped_build:
    input:
        auspice_json=auspice_dir
        + f"/{{build_name}}/latest/{auspice_prefix}_{{build_name}}.json",
        tip_freq_json=auspice_dir
        + f"/{{build_name}}/latest/{auspice_prefix}_{{build_name}}_tip-frequencies.json",
        root_sequence_json=auspice_dir
        + f"/{{build_name}}/latest/{auspice_prefix}_{{build_name}}_root-sequence.json",
    output:
        auspice_json=auspice_dir
        + f"/{{build_name}}/{{date}}/{auspice_prefix}_{{build_name}}_{{date}}.json",
        tip_freq_json=auspice_dir
        + f"/{{build_name}}/{{date}}/{auspice_prefix}_{{build_name}}_{{date}}_tip-frequencies.json",
        root_sequence_json=auspice_dir
        + f"/{{build_name}}/{{date}}/{auspice_prefix}_{{build_name}}_{{date}}_root-sequence.json",
    params:
        in_dir=auspice_dir + f"/{{build_name}}/latest/*",
        out_dir=auspice_dir + f"/{{build_name}}/{{date}}",
    wildcard_constraints:
        date="[-\d]+",
    shell:
        """
        cp {input.auspice_json} {output.auspice_json}; \
        cp {input.tip_freq_json} {output.tip_freq_json}; \
        cp {input.root_sequence_json} {output.root_sequence_json};
        """
