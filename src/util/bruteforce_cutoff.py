import os
import math
import sys

import numpy as np
from pathlib import Path
from subprocess import run, DEVNULL, TimeoutExpired

from matrix_loader import load
from recommend_cutoff import recommend


def count_pairs(path):
    if not path.exists():
        return 0

    cnt = 0
    with open(path) as fin:
        fin_iterator = iter(fin)
        next(fin_iterator)
        for line in fin_iterator:
            if ',' in line:
                cnt += 1
            elif line[0] == '#':
                break
    return cnt


delta = 1  # difference between thresholds
step = 0.5


def bruteforce(path, fasta_path=None):
    if not path.exists():
        print("Path", path, "does not exist")
        sys.exit(1)

    if fasta_path is None:
        fasta_path = path.with_suffix('.fasta')

    data = load(str(path))
    minmax = (data.ravel().min(), data.ravel().max())
    print("min={}, max={}".format(*minmax))

    recommended_cutoffs = recommend(data)
    del data

    temp_dir = Path.cwd() / "temp"
    temp_dir.mkdir(exist_ok=True)

    rounded_min = np.ceil(minmax[0] / step) * step
    binding_cutoffs = np.r_[
        np.arange(rounded_min, minmax[1] - delta, step=step), recommended_cutoffs[0]]
    nonbinding_cutoffs = np.r_[
        (binding_cutoffs + delta)[:-1], recommended_cutoffs[1]]

    best_cutoff_idx = 0
    max_num_pairs = -1

    num_pairs_list = []

    for idx, (binding_cutoff, nonbinding_cutoff) in enumerate(zip(binding_cutoffs, nonbinding_cutoffs)):
        output_path = (
            temp_dir / "cutoff_bruteforce_{}.txt".format(idx)).relative_to(Path.cwd())
        try:
            run([r'..\..\build-win\solver', str(path), '--binding-cutoff={}'.format(
                binding_cutoff), '--nonbinding-cutoff={}'.format(nonbinding_cutoff), '--out-name={}'.format(output_path), '--fasta-name={}'.format(fasta_path)], stdout=DEVNULL, stderr=DEVNULL, timeout=5)
        except TimeoutExpired:
            print('Timeout expired for cutoffs', binding_cutoff,
                  nonbinding_cutoff, file=sys.stderr)

        n_pairs = count_pairs(output_path)
        print('Cutoffs', binding_cutoff, nonbinding_cutoff,
              ':', n_pairs, 'pairs', file=sys.stderr)

        if n_pairs > max_num_pairs:
            max_num_pairs = n_pairs
            best_cutoff_idx = idx
        elif n_pairs != 0:
            num_pairs_list.append(n_pairs)
            if len(num_pairs_list) >= 3:
                if num_pairs_list[-1] < num_pairs_list[-2] < num_pairs_list[-3]:
                    print('Last 3 scores were decreasing, aborting search',
                          file=sys.stderr)
                    break

    best_cutoff = binding_cutoffs[best_cutoff_idx]

    return (best_cutoff, best_cutoff + delta)


if __name__ == "__main__":
    path = Path(sys.argv[1]) if len(
        sys.argv) > 1 else Path(r'..\..\data\PNIC.bin')

    cutoffs = bruteforce(path)
    print(*cutoffs)
