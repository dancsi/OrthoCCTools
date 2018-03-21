from sys import argv, exit
from matplotlib import pyplot as plt
import numpy as np
from matrix_loader import load


def print_help_and_exit():
    print("""USAGE: compare_matrices.py MAT1 MAT2""")
    exit(1)


if len(argv) < 3:
    print_help_and_exit()

m1 = load(argv[1])
m2 = load(argv[2])

d = np.abs(m1 - m2)

print(d.max())
plt.imshow(d, cmap='hot')
plt.show()
