import argparse
import itertools

parser = argparse.ArgumentParser(description="Generate a library of candidate peptides.")
parser.add_argument('-n', type=int, metavar='num', default=4, help="Number of heptads.")

args = parser.parse_args()
num_heptads = args.n

heptads = [''.join(['A', g, a, 'A', 'A', 'L', e]) for g in 'EK' for a in 'IN' for e in 'EK']
peptides = itertools.product(heptads, repeat=num_heptads)

for i,p in enumerate(peptides, start=1):
    seq = ''.join(p)
    print(f">P{i}\n{seq}")
