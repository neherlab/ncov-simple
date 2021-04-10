## Simplified nextstrain/ncov (SARS-CoV-2) workflow using reference data sets

### Set up

```
git clone https://github.com/nextstrain/ncov-simple.git
./decompress_example.sh
```

### make your own build
Copy the `example_profile` folder and its content and modify the file `builds.yaml`.
```
builds:
  us-washington:
    reference_metadata: "s3://nextstrain-data/ncov-intermediates/us-washington_metadata.tsv.xz"
    reference_sequences: "s3://nextstrain-data/ncov-intermediates/us-washington_alignment.fasta.xz"
    user_metadata: "example_data/metadata.tsv"
    user_sequences: "example_data/sequences.fasta"
```
You want to change the entries `user_metadata` and `user_sequences` to point at your data. 
These entries can be a list of multiple files.
Similarly, you should choose an appropriate reference data set. 

