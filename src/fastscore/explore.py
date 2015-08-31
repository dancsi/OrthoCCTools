import os, math

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm

weak_percentile = 0.21450289100201259
strong_percentile = 0.0033595580728521786

path = 'C:\\Dev\\temp\\penta.bin'
sz = os.path.getsize(path)
n_peptides = int(math.sqrt(sz / 4))
data = np.memmap(path, dtype='float32', mode='r', shape=(n_peptides, n_peptides))

def hist(data, nbins=100):
    data = data.flatten()
    h,bins=np.histogram(data, bins=nbins, density=True)
    minmax = (min(bins), max(bins))
    plt.bar(bins[:-1], h, width=(minmax[1] - minmax[0])/nbins)
    plt.xlim(minmax)

    mu, std = norm.fit(data)

    x = (bins[1:] + bins[:-1])/2
    gauss = norm.pdf(x, mu, std)
    plt.plot(x, gauss, 'r')
    
    plt.show()

    weak_cutoff = norm.ppf(weak_percentile, mu, std)
    strong_cutoff = norm.ppf(strong_percentile, mu, std)
    
    return weak_cutoff, strong_cutoff

cutoffs = hist(data)
print(cutoffs)
