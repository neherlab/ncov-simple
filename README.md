**This repository currently serves as a sandbox to explore how our main [ncov workflow](github.com/nextstrain/ncov) could restructured and simplified. Feel free to use it, but don't expect continuity or extensive documentation at this point.**

## Simplified nextstrain/ncov (SARS-CoV-2) workflow using reference data sets

The reference data sets are assembled from SARS-CoV-2 available in the public domain.
The availability of open data varies greatly be region and the only regions for which we have good coverage at the moment are the United States and Oceania.

### Set up

```
git clone https://github.com/nextstrain/ncov-simple.git
./decompress_example.sh
```
The second step is only necessary if you want to run the example profile.
The standard nextstrain conda environment should be sufficient to run this workflow.
It can then be run as
```
snakemake --profile example_profile
```

### make your own build
Copy the `example_profile` folder and its content and modify the file `builds.yaml`.
The example profile contains builds that are specified like this example for an extra small build using a reference set for Washington, US:
```
builds:
  us-washington:
    reference_metadata:  "https://data.nextstrain.org/ncov-reference-datasets/US-WA-xsmall_metadata.tsv.xz"
    reference_sequences: "https://data.nextstrain.org/ncov-reference-datasets/US-WA-xsmall_sequences.fasta.xz"
    user_metadata: "example_data/metadata.tsv"
    user_sequences: "example_data/sequences.fasta"
```
You want to change the entries `user_metadata` and `user_sequences` to point at your data and choose and appropriate reference data set.
These entries can be a list of multiple files.


### Available reference data sets

Currently, we provide reference datasets representing all available data, the standard Nextstrain regions, and the US states.
Each data set comes in three sizes -- one standard size with 4k sequences, a small one with 2k sequences, and an extra small one with 1k sequences.
The sequences are a random sample from the focal region (emphasis on the last 4 months) and augmented by sequences from outside the focal regions using Nextstrain's proximity guided sampling.

The data are available at `data.nextstrain.org/ncov-reference-datasets` and follow the URL patter

```
https://data.nextstrain.org/ncov-reference-datasets/<region>-<size>_metadata.tsv.xz"
https://data.nextstrain.org/ncov-reference-datasets/<region>-<size>_sequences.fasta.xz"
```
where `region` can be `global`, `africa`, `asia`, `europe`, `north-america`, `oceania`, `south-america` or a US state using two letter state abbreviations as in `US-CA`, `US-NJ`, etc.


## Templated builds

In some situation, you want to run builds for many different data subsamples, for example for different countries, states, variants etc.
To facilitate such builds, this repo contains a functionality to __template__ builds.
This is used in the profiles `profiles/basel-swiss`, `profiles/basel-countries`, and `profiles/genbank` and work as follows:

 * you specify lists of features for which you want to make builds. For example regions, countries, build sizes
 * you specify a number of parameters for the builds, these are expressions that will be evaluated
 * you specify metadata adjustments, which are pandas queries
 * you specify subsampling schemes.

The builds spec would then look like this:
```
templated-builds:
  regions:
    build_patterns:
      region:
        north-america: "North America"
        asia: "Asia"
        africa: "Africa"
        europe: "Europe"
        south-america: "South America"
        oceania: "Oceania"
      size:
        4k: 4000
        2k: 2000
        1k: 1000
    build_name: "{region}_{size}"
    ...
```
This defines a templated build `regions` that will generate builds of three different sizes for each region.
The subsampling parameters are defined as
```
    subsampling_parameters:
      date_cutoff: "(datetime.date.today() - datetime.timedelta(weeks=18)).strftime('%Y-%m-%d')"
      s_global_early: "int(0.1*{size})"
      s_global_late:  "int(0.2*{size})"
      s_region_early: "int(0.3*{size})"
      s_region_late:  "int(0.4*{size})"
    subsamples:
      global_early:
        filters: "--exclude-where region='{region}' --group-by country year month --subsample-max-sequences {s_global_early} --max-date {date_cutoff}"
      region_early:
        filters: "--exclude-where region!='{region}' --group-by country year month --subsample-max-sequences {s_region_early} --max-date {date_cutoff}"
      global_late:
        filters: "--exclude-where region='{region}' --group-by country year month --subsample-max-sequences {s_global_late} --min-date {date_cutoff}"
        priorities: "region_late"
      region_late:
        filters: "--exclude-where region!='{region}' --group-by country year month --subsample-max-sequences {s_region_late} --min-date {date_cutoff}"
```
Each in each expression, the variables are substituted before the command is launched.
Metadata adjustments (for example to replace all divisions outside a focal region with a higher order geographic designation) can be implemented as
```
    metadata_adjustments:
      - query: region!='{region}'
        dst: country
        src: region
      - query: region!='{region}'
        dst: division
        src: region
      - query: region!='{region}'
        dst: location
        src: region
```
Again, the build parameters and subsampling parameters will be substituted before the expressions are evaluated.


