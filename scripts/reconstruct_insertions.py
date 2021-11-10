#! /usr/local/Caskroom/mambaforge/base/envs/nextstrain/bin/python

# Get all insertions and code as binary vector: Use metadata for this
# Load tree
# Reconstruct using tree time

import json

import augur.ancestral as ancestral
import Bio.Align
import click
import pandas as pd


@click.command()
@click.option("--metadata", default="builds-test/test/metadata.tsv", type=str)
@click.option("--tree", default="builds-test/test/tree.nwk", type=str)
@click.option("--output", default="builds-test/test/insertions.json", type=click.File("w"))
def main(metadata, tree, output):
    meta = pd.read_csv(metadata, sep="\t", index_col=0)
    # Find all unique insertions, create list
    # While waiting for insertions, use pango for fun
    # Tricky: Do separately or jointly for insertions and frameshifts? Probably easier.

    def split_csv_field(x):
        """Split csv field into list, with empty strings yielding []"""
        return x.split(",") if type(x) == str else []

    # Construct of insertions
    characters = {
        insertion for insertions in meta.insertions for insertion in split_csv_field(insertions)
    }

    # Dictionary that maps property to position in pseudosequence
    mapping = {}
    for pos, character in enumerate(characters):
        mapping[character] = pos

    inverse_mapping = {v: k for k, v in mapping.items()}

    # Take list of characters and turn into binary vector
    # Convention: A = 0, G = 1
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
    print(meta.insertions_vector[meta.insertions_vector.apply(lambda x: "G" in x) > 0])

    # Transform binary vector into BioSeq alignment
    alignment = Bio.Align.MultipleSeqAlignment(
        [
            Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(vector), id=strain, name=strain)
            for strain, vector in meta.insertions_vector.iteritems()
        ]
    )
    print(alignment)

    from treetime import TreeAnc

    tt = TreeAnc(tree=tree, aln=alignment)
    tt.infer_ancestral_sequences()
    nodes_with_mutations = ancestral.collect_mutations_and_sequences(tt)

    def mut_to_str(mut: str) -> str:
        """Convert mutation to string"""
        indel_type = "ins" if mut[0] == "A" else "rev ins"
        return f"{indel_type} {inverse_mapping[int(mut[1:-1])-1]}"

        return mut.replace("-", "")

    for node in nodes_with_mutations["nodes"]:
        nodes_with_mutations["nodes"][node]["muts"] = list(
            map(mut_to_str, nodes_with_mutations["nodes"][node]["muts"])
        )
        if nodes_with_mutations["nodes"][node]["muts"] != []:
            print(f"{node}: {nodes_with_mutations['nodes'][node]['muts']}")

    json.dump(nodes_with_mutations, output)


if __name__ == "__main__":
    main()
