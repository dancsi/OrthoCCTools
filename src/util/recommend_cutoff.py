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

path = Path(sys.argv[1] if len(sys.argv) >
            1 else r"""..\..\data\dipeptide-at-least-four.bin""")

if not path.exists():
    print("Path", path, "does not exist")
    sys.exit(1)

data = load(str(path))


def hist(data, nbins=100):
    data = data.flatten()
    h, bins = np.histogram(data, bins=nbins, density=True)
    minmax = (min(bins), max(bins))
    plt.bar(bins[:-1], h, width=(minmax[1] - minmax[0]) / nbins)
    plt.xlim(minmax)

    mu, std = norm.fit(data)
    real_dist = rv_histogram((h, bins))

    x = (bins[1:] + bins[:-1]) / 2
    gauss = norm.pdf(x, mu, std)
    plt.plot(x, gauss, 'k.')

    weak_cutoff = real_dist.ppf(weak_percentile)
    strong_cutoff = real_dist.ppf(strong_percentile)

    weak_cutoff_gaussian = norm.ppf(weak_percentile, mu, std)
    strong_cutoff_gaussian = norm.ppf(strong_percentile, mu, std)

    plt.axvline(x=weak_cutoff, color='r')
    plt.axvline(x=strong_cutoff, color='r')

    plt.axvline(x=weak_cutoff_gaussian, color='k')
    plt.axvline(x=strong_cutoff_gaussian, color='k')

    plt.show()

    return weak_cutoff, strong_cutoff


cutoffs = hist(data)
print(cutoffs)
