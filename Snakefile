import datetime


localrules:
    all,
    clean,
    clean_all,
    deploy_nodate,
    deploy_date,
    deploy_all,
    dump_config,  #, deploy_force, deploy_all_force


if "builds" not in config:
    config["builds"] = {}

if "files" not in config:

    configfile: "defaults/parameters.yaml"


if "origins" in config:

    include: "workflow/snakemake_rules/preprocess.smk"


if "reference-builds" in config:
    config["builds"].update(config["reference-builds"])

    # Include rules to handle primary build logic from multiple sequence alignment
    # to output of auspice JSONs for a default build.
    include: "workflow/snakemake_rules/reference_build.smk"


if "templated-builds" in config:

    include: "workflow/snakemake_rules/templated_build.smk"


if len(config["builds"]):

    include: "workflow/snakemake_rules/subsampling.smk"
    include: "workflow/snakemake_rules/core.smk"


auspice_prefix = config.get("auspice_prefix", "ncov")
auspice_dir = config.get("auspice_dir", "auspice")
build_dir = config.get("build_dir", "builds")

date = datetime.date.today()
suffixes = ["", "_root-sequence", "_tip-frequencies"]


rule all:
    input:
        [
            f"{auspice_dir}/{build}/latest/{auspice_prefix}_{build}{suffix}.json"
            for build in config["builds"]
            for suffix in suffixes
        ],
    params:
        slack_hook=config.get("slackHook", "google.com"),
    shell:
        """
        curl -X POST -H 'Content-type: application/json' \
        --data '{{"text":"Builds done, ready for deployment"}}' \
        {params.slack_hook}
        """


# def deploy_files(w):
#     return " ".join([f"{auspice_dir}/{w.build}/{auspice_prefix}_{w.build}{w.date}{suffix}.json" for suffix in suffixes])

# rule deploy_force:
#     # input: ancient([f"{auspice_dir}/{auspice_prefix}_{{build}}{{date}}{suffix}.json" for suffix in suffixes])
#     output: 'deploy-force/{build,[^_]+}{date,.{0}|_.+}_force',
#     # nexde url1 input; nexde url2 input
#     params: lambda w, input, output: " ; ".join([f'nextstrain deploy {url} {deploy_files(w)} 2>&1 | tee -a {output}' for url in config["builds"][w.build]["deploy_urls"]])
#     shell: '{params}'


rule deploy_nodate:
    input:
        [
            f"{auspice_dir}/{{build}}/latest/{auspice_prefix}_{{build}}{suffix}.json"
            for suffix in suffixes
        ],
    output:
        "deploy/{build}/latest",
    # nexde url1 input; nexde url2 input
    params:
        lambda w, input, output: " ; ".join(
            [
                f"nextstrain remote upload {url}/{w.build} {input} 2>&1 | tee -a {output}"
                for url in config["builds"][w.build]["deploy_urls"]
            ]
        ),
    shell:
        "{params}"


rule deploy_date:
    input:
        [
            auspice_dir
            + f"/{{build}}/{{date}}/{auspice_prefix}_{{build}}_{{date}}{suffix}.json"
            for suffix in suffixes
        ],
    output:
        "deploy/{build}/{date,[\d-]+}",
    # nexde url1 input; nexde url2 input
    params:
        lambda w, input, output: " ; ".join(
            [
                f"nextstrain remote upload {url}/{w.build}_{w.date} {input} 2>&1 | tee -a {output}"
                for url in config["builds"][w.build]["deploy_urls"]
            ]
        ),
    shell:
        "{params}"


rule deploy_all:
    input:
        expand("deploy/{build}/latest", build=config["builds"]),
        expand("deploy/{build}/{date}", build=config["builds"], date=date),


# rule deploy_all_force:
#     input:
#         expand("deploy/{build}_force", build=config["builds"])


def dated_deploy(build_name: str):
    return f"deploy/{build_name}/{date}"


def undated_deploy(build_name: str):
    return f"deploy/{build_name}/latest"


deploy_test_list = [
    "europe",
    "switzerland",
    "germany",
    "na-usa",
    "united-kingdom",
    "global",
    "south-america",
    "north-america",
    "asia",
    "africa",
    "oceania",
]


rule deploy_test:
    input:
        [dated_deploy(build_name) for build_name in deploy_test_list],
        [undated_deploy(build_name) for build_name in deploy_test_list],


rule clean_all:
    params:
        "pre-processed",
        "data",
        "log/*",
        "logs/*",
        "benchmarks",
        auspice_dir,
        build_dir,
        "deploy/*",
        "freezed",
        "builds*",
        "auspice*",
    shell:
        "rm -rfv {params}"


rule clean:
    params:
        "log/*",
        "logs/*",
        "benchmarks",
        auspice_dir,
        build_dir,
        "deploy/*",
        "freezed",
    shell:
        "rm -rfv {params}"


rule dump_config:
    run:
        import sys
        import ruamel.yaml

        yaml = ruamel.yaml.YAML()
        yaml.dump(config["builds"], sys.stdout)
