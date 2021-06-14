import numpy as np
from treetime import TreeAnc
from Bio import AlignIO
from collections import defaultdict

comparison_seqs = {s.id: np.array(s) for s in AlignIO.read('../2021-05_origins-letter/rooting_and_trees/sarbeco.aligned.fasta', 'fasta')}

tt = TreeAnc(tree = "builds-early/early-19B/tree_raw.nwk", aln="builds-early/early-19B/masked.fasta", gtr='JC69')

tt.infer_ancestral_sequences(reconstruct_tip_states=True)

distance = defaultdict(dict)
comp_seq_name = 'hCoV-19/bat/Yunnan/RaTG13/2013|EPI_ISL_402131|2013-07-24'
short_name = comp_seq_name.split('|')[0].split('/')[-2]
comp_seq = comparison_seqs[comp_seq_name]
mindist = 30000
mindist_terminal = 30000
for n in tt.tree.find_clades():
    dist = np.sum(tt.sequence(n, reconstructed=True, as_string=False)!=comp_seq)
    if dist<=mindist:
        print(n.name, dist)
        mindist=dist
    if dist<=mindist_terminal:
        print(n.name, dist)
        mindist_terminal=dist
    distance[n.name][f"distance_{comp_seq_name}"] =  dist

import json

with open("bat_distances.json", 'w') as fh:
    json.dump({"nodes": distance}, fh)

key = f"distance_{comp_seq_name}"
with open("bat_distances.tsv", 'w') as fh:
    fh.write('strain\tdistance3\n')
    for n, v in sorted(distance.items(), key=lambda x:x[1][key]):
        if not n.startswith('NODE'):
            fh.write(f'{n}\t{min(mindist_terminal+10, v[key])}\n')
