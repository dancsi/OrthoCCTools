from matrix_loader import load
from matplotlib import pyplot as plt
from pathlib import Path
from subprocess import run, PIPE
import math

import sys


def plot_single(path, ax):
    if not path.exists():
        print("Path", path, "does not exist")
        sys.exit(1)

    bin_path = path.with_suffix('.bin')
    fasta_path = bin_path.with_suffix('.fasta')

    if not bin_path.exists():
        run([r'..\..\build-win\fastscore', str(fasta_path)])

    mm = load(str(bin_path))
    im = ax.imshow(mm, cmap='hot')
    return im

if __name__ == "__main__":

    if len(sys.argv) < 2:
        sys.argv.append(r"""..\..\data\tetraheptads.orthoset.bin""")

    paths = [Path(p) for p in sys.argv[1:]]

    n = len(paths)
    if n <= 3:
        dimensions = [(1, 1), (1, 2), (2, 2)]
        (nrows, ncols) = dimensions[n - 1]
    else:
        divisors = [num for num in range(
            1, round(math.sqrt(n)) + 1) if n % num == 0 and num * num <= n]
        nrows = divisors[-1]
        ncols = n // nrows

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols)

    if n == 1:
        axes.flat = [axes]
        display_axes = True
    else:
        display_axes = False

    for path, ax in zip(paths, axes.flat):
        im = plot_single(path, ax)

    if not display_axes:
        for ax in axes.flat:
            ax.axis('off')

    fig.colorbar(im, ax=list(axes.flat))

    plt.show()
