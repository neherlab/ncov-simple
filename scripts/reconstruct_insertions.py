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
        vector = ['A'] * len(characters)
        for character in character_list:
            vector[mapping[character]] = 'G'
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
            Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(vector), id=strain)
            for strain, vector in meta.insertions_vector.iteritems()
        ]
    )
    print(alignment)

    print(ancestral.ancestral_sequence_inference(tree,alignment))
    """
    Bug:
    Traceback (most recent call last):
      File "/Users/cr/code/ncov-simple/scripts/reconstruct_insertions.py", line 69, in <module>
        main()
      File "/usr/local/Caskroom/mambaforge/base/envs/nextstrain/lib/python3.8/site-packages/click/core.py", line 1128, in __call__
        return self.main(*args, **kwargs)
      File "/usr/local/Caskroom/mambaforge/base/envs/nextstrain/lib/python3.8/site-packages/click/core.py", line 1053, in main
        rv = self.invoke(ctx)
      File "/usr/local/Caskroom/mambaforge/base/envs/nextstrain/lib/python3.8/site-packages/click/core.py", line 1395, in invoke
        return ctx.invoke(self.callback, **ctx.params)
      File "/usr/local/Caskroom/mambaforge/base/envs/nextstrain/lib/python3.8/site-packages/click/core.py", line 754, in invoke
        return __callback(*args, **kwargs)
      File "/Users/cr/code/ncov-simple/scripts/reconstruct_insertions.py", line 65, in main
        print(ancestral.ancestral_sequence_inference(tree,alignment))
      File "/usr/local/Caskroom/mambaforge/base/envs/nextstrain/lib/python3.8/site-packages/augur/ancestral.py", line 48, in ancestral_sequence_inference
        tt = TreeAnc(tree=tree, aln=aln, ref=ref, gtr='JC69',
      File "/usr/local/Caskroom/mambaforge/base/envs/nextstrain/lib/python3.8/site-packages/treetime/treeanc.py", line 162, in __init__
        self._check_alignment_tree_gtr_consistency()
      File "/usr/local/Caskroom/mambaforge/base/envs/nextstrain/lib/python3.8/site-packages/treetime/treeanc.py", line 380, in _check_alignment_tree_gtr_consistency
        raise MissingDataError("TreeAnc._check_alignment_tree_gtr_consistency: At least 30\\% terminal nodes cannot be assigned a sequence!\n"
    treetime.MissingDataError: TreeAnc._check_alignment_tree_gtr_consistency: At least 30\% terminal nodes ca
    """


if __name__ == "__main__":
    main()
