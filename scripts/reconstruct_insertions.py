#! /usr/local/Caskroom/mambaforge/base/envs/nextstrain/bin/python

# Get all insertions and code as binary vector: Use metadata for this
# Load tree
# Reconstruct using tree time

import click
import pandas as pd


@click.command()
@click.option("--metadata", default="builds-test/test/metadata.tsv", type=click.File("r"))
def main(metadata):
    meta = pd.read_csv(metadata, sep="\t", index_col=0)
    # Find all unique insertions, create list
    # While waiting for insertions, use pango for fun
    # Tricky: Do separately or jointly for insertions and frameshifts? Probably easier.

    def split_csv_field(x):
        """Split csv field into list, with empty strings yielding []"""
        return x.split(",") if type(x) == str else []

    characters = {
        insertion for insertions in meta.insertions for insertion in split_csv_field(insertions)
    }

    # Dictionary that maps property to position in pseudosequence
    mapping = {}
    for pos, character in enumerate(characters):
        mapping[character] = pos

    # Take list of characters and turn into binary vector
    def character_list_to_vector(character_list_csv):
        """Turn comma separated list into binary vector"""
        character_list = split_csv_field(character_list_csv)
        vector = [0] * len(characters)
        for character in character_list:
            vector[mapping[character]] = 1
        return vector

    meta = meta.assign(insertions_vector=meta.insertions.apply(character_list_to_vector))
    print(meta.insertions_vector[0])


if __name__ == "__main__":
    main()
