from itertools import chain
from Bio import SeqIO

fpairs = open('pairs.txt')
pairs = list(map( lambda t: (int(t[0][1:]), int(t[1][1:])), [tuple(line.strip().split(',')) for line in fpairs]))
seqs = dict.fromkeys(chain.from_iterable(pairs), None)

for record in enumerate(SeqIO.parse(open('../data/dipeptide-at-least-four.fasta'), 'fasta')):
    if record[0] in seqs:
        seqs[record[0]] = str(record[1].seq)

fout = open('transformed.txt', "w")

for pair in pairs:
    print("%s,%s"%(seqs[pair[0]], seqs[pair[1]]), file=fout)

fout.close()
