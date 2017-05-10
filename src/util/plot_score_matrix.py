from matrix_loader import load
from matplotlib import pyplot as plt
from pathlib import Path
from subprocess import run, PIPE

import sys

path = Path(sys.argv[1] if len(sys.argv) >
            1 else r'..\..\data\dipeptide-at-least-four.orthoset.fasta')

bin_path = path.with_suffix('.bin')
fasta_path = bin_path.with_suffix('.fasta')

if not bin_path.exists():
    run([r'..\..\build-win\fastscore', str(fasta_path)])

mm = load(str(bin_path))

fig, ax1 = plt.subplots()
cax = ax1.imshow(mm, cmap='hot')
fig.colorbar(cax)

# ax2.hist(mm.ravel())

plt.show()
