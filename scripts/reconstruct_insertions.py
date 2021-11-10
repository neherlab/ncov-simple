#! /usr/local/Caskroom/mambaforge/base/envs/nextstrain/bin/python

# Get all insertions and code as binary vector: Use metadata for this
# Load tree
# Reconstruct using tree time

import augur.ancestral as ancestral
import Bio.Align
import click
import pandas as pd


@click.command()
@click.option("--metadata", default="builds-test/test/metadata.tsv", type=str)
@click.option("--tree", default="builds-test/test/tree.nwk", type=str)
def main(metadata, tree):
    meta = pd.read_csv(metadata, sep="\t", index_col=0)
    # Find all unique insertions, create list
    # While waiting for insertions, use pango for fun
    # Tricky: Do separately or jointly for insertions and frameshifts? Probably easier.

    def split_csv_field(x):
        """Split csv field into list, with empty strings yielding []"""
        return x.split(",") if type(x) == str else []

    print(meta.columns)
    # Construct of insertions
    characters = {
        insertion for insertions in meta.insertions for insertion in split_csv_field(insertions)
    }

    # Dictionary that maps property to position in pseudosequence
    mapping = {}
    for pos, character in enumerate(characters):
        mapping[character] = pos

    # Take list of characters and turn into binary vector
    # Convention: A = 0, C = 1
    def character_list_to_vector(character_list_csv):
        """Turn comma separated list into binary vector"""
        character_list = split_csv_field(character_list_csv)
        vector = ["A"] * len(characters)
        for character in character_list:
            vector[mapping[character]] = "G"
        return "".join(vector)

    meta = meta.assign(insertions_vector=meta.insertions.apply(character_list_to_vector))

    meta.to_csv("builds-test/test/metadata_with_insertions.tsv", sep="\t")

    # Print mapping
    print(mapping)
    # Print all sequences with insertions
    print(meta.insertions_vector[meta.insertions_vector.apply(lambda x: "C" in x) > 0])

    # Transform binary vector into BioSeq alignment
    alignment = Bio.Align.MultipleSeqAlignment(
        [
            Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(vector), id=strain, name=strain)
            for strain, vector in meta.insertions_vector.iteritems()
        ]
    )

    from treetime import TreeAnc

    tt = TreeAnc(tree=tree, aln=alignment)
    tt.infer_ancestral_sequences()


if __name__ == "__main__":
    main()
