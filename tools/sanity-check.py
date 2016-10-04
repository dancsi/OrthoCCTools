import argparse
import csv
import os
import subprocess
from itertools import chain

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("orthoset",
                    type=str,
                    nargs="?",
                    default="pairs.txt",
                    help="The file containing comma-separated peptide ids of the pairs, output by the solver")

parser.add_argument("fasta",
                    type=str,
                    nargs="?",
                    default="../data/dipeptide-at-least-four.fasta",
                    help="The file containing the peptides")

args = parser.parse_args()

pairs = set()
with open(args.orthoset, 'r', newline='') as fpairs:
    reader = csv.reader(fpairs)
    for row in reader:
        (id1, id2) = (row[0][1:], row[1][1:])
        (i, j) = (int(id1), int(id2))
        if i > j:
            i, j = j, i
        pairs.add((i, j))

seqs = dict.fromkeys(chain.from_iterable(pairs), None)

with open('temp.fasta', 'w') as fout:
    for record in enumerate(SeqIO.parse(open(args.fasta), 'fasta')):
        if record[0] in seqs:
            print(">{}\n{}".format(record[0], record[1].seq), file=fout)

output = subprocess.run(['python', 'bzipscore-all.py', 'temp.fasta'],
                        stdout=subprocess.PIPE) \
    .stdout.decode('utf-8').strip().splitlines()

min_binding = float('inf')
max_nonbinding = 0

reader = csv.reader(output)
for line in reader:
    (id1, id2, score) = tuple(line)
    (i, j, s) = (int(id1), int(id2), float(score))
    if i > j:
        i, j = j, i
    if (i, j) in pairs:
        if abs(s) < abs(min_binding):
            min_binding = s
    else:
        if abs(s) > abs(max_nonbinding):
            max_nonbinding = s

print("min binding:", min_binding)
print("max nonbinding:", max_nonbinding)

os.remove('temp.fasta')
