import argparse
import csv
import math
import os
import sys
from collections import namedtuple
from pathlib import Path
from subprocess import DEVNULL, TimeoutExpired, run

import numpy as np

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


BruteforceResult = namedtuple("BruteforceResult", [
                              "binding_cutoff", "nonbinding_cutoff", "binding_cutoffs", "nonbinding_cutoffs", "num_pairs", "best_cutoff_idx", "temp_dir"])


def bruteforce(path, fasta_path=None, fast_exit=True, delta=1, step=0.5, timeout=5):
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
                binding_cutoff), '--nonbinding-cutoff={}'.format(nonbinding_cutoff), '--out-name={}'.format(output_path), '--fasta-name={}'.format(fasta_path)], stdout=DEVNULL, stderr=DEVNULL, timeout=timeout)
        except TimeoutExpired:
            print('Timeout expired for cutoffs', binding_cutoff,
                  nonbinding_cutoff, file=sys.stderr)

        n_pairs = count_pairs(output_path)
        print('Cutoffs', binding_cutoff, nonbinding_cutoff,
              ':', n_pairs, 'pairs', file=sys.stderr)

        if n_pairs > max_num_pairs:
            max_num_pairs = n_pairs
            best_cutoff_idx = idx
        if n_pairs != 0:
            num_pairs_list.append(n_pairs)
            if fast_exit and len(num_pairs_list) >= 3:
                if num_pairs_list[-1] < num_pairs_list[-2] < num_pairs_list[-3]:
                    print('Last 3 scores were decreasing, aborting search',
                          file=sys.stderr)
                    break

    return BruteforceResult(binding_cutoff=binding_cutoffs[best_cutoff_idx],                            nonbinding_cutoff=nonbinding_cutoffs[best_cutoff_idx],
                            binding_cutoffs=binding_cutoffs,
                            nonbinding_cutoffs=nonbinding_cutoffs,
                            num_pairs=num_pairs_list,
                            best_cutoff_idx=best_cutoff_idx,
                            temp_dir=temp_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('path', type=Path)
    parser.add_argument('--delta', default=1, type=float)
    parser.add_argument('--step', default=0.5, type=float)
    parser.add_argument('--fast_exit', action='store_true')
    parser.add_argument('--timeout', default=5, type=float)
    args = parser.parse_args()

    cutoffs = bruteforce(args.path, delta=args.delta, step=args.step,
                         fast_exit=args.fast_exit, timeout=args.timeout)

    print(cutoffs.nonbinding_cutoff, cutoffs.binding_cutoff)
    print(cutoffs.nonbinding_cutoffs)
    print(cutoffs.binding_cutoffs)
    print(cutoffs.num_pairs)

    csv_rows = zip(cutoffs.nonbinding_cutoffs, cutoffs.binding_cutoffs,
                   cutoffs.num_pairs)

    with open(cutoffs.temp_dir / "cutoff_bruteforce.csv", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerows(csv_rows)
