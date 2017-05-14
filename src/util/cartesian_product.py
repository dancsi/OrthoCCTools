import sys
from Bio import SeqIO
from pathlib import Path
from itertools import product


def product_fasta(set_paths, dest_path):
    sets = list(
        map(lambda t: (list(SeqIO.parse(str(t[0]), 'fasta')), t[1]), set_paths))

    current_length = 0

    factors = []

    for i, (sequences, cnt) in enumerate(sets):
        for rep in range(cnt):
            characters_to_cut = current_length % 7

            def transformer(x): return x
            if characters_to_cut != 0:
                if characters_to_cut != 1:
                    print(
                        'Currently we only support sequences that begin with at most one "-"')
                    sys.exit(1)

                def transformer(x): return x[1:]

            factor = [transformer(str(record.seq)) for record in sequences]
            current_length += len(transformer(str(sequences[0].seq)))

            factors.append(factor)

    prod = product(*factors)
    new_sequences = map(lambda t: ''.join(t), prod)

    with open(dest_path, 'w') as fout:
        cnt = 0
        for seq in new_sequences:
            cnt += 1
            print('>P{}\n{}'.format(cnt, seq), file=fout)
        print(cnt, 'sequences generated')


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(r'''Construct the cartesian product of fasta files''')
        print(
            r'''USAGE: cartesian_product.py dest_path path1 [-n NUM] [path2 [-n NUM] ...]''')
        sys.exit(1)

    dest_path = Path(sys.argv[1])
    n_args = len(sys.argv)
    i = 2
    set_paths = []

    while i < n_args:
        path = Path(sys.argv[i])
        if not path.exists():
            print('Path', path, 'does not exist')
            sys.exit(1)
        num = 1
        if i + 1 < n_args and sys.argv[i + 1] == '-n':
            num = int(sys.argv[i + 2])
            i += 2
        set_paths.append((path, num))
        i += 1

    product_fasta(set_paths, dest_path)
