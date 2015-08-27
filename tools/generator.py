import itertools, argparse, random, sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Generate some peptides")
parser.add_argument('-n', type=int, metavar='num', default=None, help="Number of peptides. If not present, all possible will be generated")
parser.add_argument('-np', type=int, metavar='num', default=4, help="Number of hepdates.")

#parser.add_argument('--full', action='store_true', default=False, help="Use the full interaction set (off by default)")
args = parser.parse_args()

n = args.n;
np = args.np;

def iter_sample_fast(iterable, samplesize):
    results = []
    iterator = iter(iterable)
    # Fill in the first samplesize elements:
    for _ in range(samplesize):
        results.append(next(iterator))
    random.shuffle(results)  # Randomize their positions
    for i, v in enumerate(iterator, samplesize):
        r = random.randint(0, i)
        if r < samplesize:
            results[r] = v  # at a decreasing rate, replace random items

    if len(results) < samplesize:
        raise ValueError("Sample larger than population.")
    return results

heptads = [''.join(['A', g, a, 'A', 'A', 'L', e]) for g in 'EK' for a in 'IN' for e in 'EK']

peptides = itertools.product(heptads, repeat=np)
if n!=None:
    peptides = random.sample(list(peptides), n)
    #peptides = iter_sample_fast(peptides, n)
peptides = map(lambda p: ''.join(p), peptides)
records = map(lambda p: SeqRecord(Seq(p[1], generic_protein), id="P"+str(p[0]), description=""), enumerate(peptides, start=1))

SeqIO.write(records, sys.stdout, "fasta")
