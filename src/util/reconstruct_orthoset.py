from itertools import chain
from Bio import SeqIO
import sys
from pathlib import Path


pairs_path = Path(sys.argv[1] if len(sys.argv) >
                  1 else r'..\..\data\dipeptide-at-least-four.pairs')
fasta_path = pairs_path.with_suffix('.fasta')
orthoset_fasta_path = pairs_path.with_suffix('.orthoset.fasta')

peptides = []

with open(pairs_path) as fin:
    fin_iterator = iter(fin)
    next(fin_iterator)
    for line in fin_iterator:
        if ',' in line:
            p1, p2 = line.strip().split(',')
            peptides.append(p1)
            peptides.append(p2)
        elif line[0] == '#':
            break

peptide_set = set(peptides)
sequences = dict()

for record in SeqIO.parse(str(fasta_path), 'fasta'):
    if record.id in peptide_set:
        sequences[record.id] = record.seq

with open(orthoset_fasta_path, 'w') as fout:
    prev = None
    for peptide_id in peptides:
        if peptide_id != prev:
            seq = sequences[peptide_id]
            print('>{}\n{}'.format(peptide_id, seq), file=fout)
            prev = peptide_id
