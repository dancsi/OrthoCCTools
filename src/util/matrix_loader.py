import numpy as np
from matplotlib import pyplot as plt


def load(path):
    header = np.memmap(path, dtype=np.uint64, mode='r')
    n, m = header[1], header[2]
    offset = int(header[3])
    print("Got offset", offset)
    del header
    return np.memmap(path, dtype='float32', mode='r', offset=offset, shape=(n, m))


mm = load("../../build-win/chained.bin")

fig, ax = plt.subplots()
cax = plt.imshow(mm, cmap='hot')
fig.colorbar(cax)
plt.show()
