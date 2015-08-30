import numpy as np
from matplotlib import pyplot as plt

path = 'C:\\Dev\\temp\\penta.bin'
data = np.memmap(path, dtype='float32', mode='r', shape=(1<<15, 1<<15))


