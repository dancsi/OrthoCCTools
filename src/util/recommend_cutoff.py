import os
import math
import sys

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm, rv_histogram
from matrix_loader import load
from pathlib import Path

weak_percentile = 0.21450289100201259
strong_percentile = 0.0033595580728521786


def recommend(data, nbins=100, plot=False):
    data = data.flatten()
    h, bins = np.histogram(data, bins=nbins, density=True)
    minmax = (min(bins), max(bins))

    real_dist = rv_histogram((h, bins))

    weak_cutoff = real_dist.ppf(weak_percentile)
    strong_cutoff = real_dist.ppf(strong_percentile)

    mu, std = norm.fit(data)
    weak_cutoff_gaussian = norm.ppf(weak_percentile, mu, std)
    strong_cutoff_gaussian = norm.ppf(strong_percentile, mu, std)

    if plot:
        plt.bar(bins[:-1], h, width=(minmax[1] - minmax[0]) / nbins)
        plt.xlim(minmax)

        x = (bins[1:] + bins[:-1]) / 2
        gauss = norm.pdf(x, mu, std)
        plt.plot(x, gauss, 'k.')

        plt.axvline(x=weak_cutoff, color='r')
        plt.axvline(x=strong_cutoff, color='r')

        plt.axvline(x=weak_cutoff_gaussian, color='k')
        plt.axvline(x=strong_cutoff_gaussian, color='k')

        plt.show()

    return strong_cutoff, weak_cutoff


if __name__ == "__main__":
    path = Path(sys.argv[1] if len(sys.argv) >
                1 else r"""..\..\data\dipeptide-at-least-four.bin""")

    if not path.exists():
        print("Path", path, "does not exist")
        sys.exit(1)

    data = load(str(path))
    cutoffs = recommend(data, plot=True)
    print(*cutoffs)
