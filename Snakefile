import datetime

if "builds" not in config:
    config["builds"] = {}
if "files" not in config:
    configfile: "defaults/parameters.yaml"

if "origins" in config:
    include: "workflow/snakemake_rules/preprocess.smk"

if "reference-builds" in config:
    rule all:
        input:
            expand("auspice/ncov_{build_name}.json", build_name=config["reference-builds"].keys())

    config["builds"].update(config["reference-builds"])
    # Include rules to handle primary build logic from multiple sequence alignment
    # to output of auspice JSONs for a default build.
    include: "workflow/snakemake_rules/reference_build.smk"

if len(config["builds"]):
    include: "workflow/snakemake_rules/subsampling.smk"
    include: "workflow/snakemake_rules/core.smk"

rule clean_all:
    message: "Removing directories: {params}"
    params:
        "builds",
        "auspice",
        "pre-processed",
        "data"
    shell:
        "rm -rfv {params}"


rule clean:
    message: "Removing directories: {params}"
    params:
        "builds",
        "auspice"
    shell:
        "rm -rfv {params}"
