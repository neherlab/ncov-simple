import datetime

if "builds" not in config:
    config["builds"] = {}

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

include: "workflow/snakemake_rules/core.smk"

rule clean:
    message: "Removing directories: {params}"
    params:
        "results",
        "auspice"
    shell:
        "rm -rfv {params}"

rule clean_all:
    message: "Removing directories: {params}"
    params:
        "results",
        "auspice",
        "downloads"
    shell:
        "rm -rfv {params}"
