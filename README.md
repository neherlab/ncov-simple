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
