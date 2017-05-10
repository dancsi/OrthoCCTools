import numpy as np
from matplotlib import pyplot as plt


def load(path):
    header = np.fromfile(path, dtype=np.uint64, count=4)
    n, m = header[1], header[2]
    offset = int(header[3])
    print("Got offset", offset)
    del header
    return np.memmap(path, dtype='float32', mode='r', offset=offset, shape=(n, m))
