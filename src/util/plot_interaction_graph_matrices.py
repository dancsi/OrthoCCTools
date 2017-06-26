from matrix_loader import load

import sys
import argparse
import networkx as nx
from pathlib import Path
from matplotlib import pyplot as plt

def plot_graph(mat, m, suffix, draw_graph=True, ax=None, edge_color='k',
               edge_width=1, edge_list=None):
    if draw_graph:
        if ax is None:
            fig, ax = plt.subplots(1, 2)
        else:
            fig = plt.gcf()
        fig.set_size_inches(11, 5)

        if edge_list is None:
            G = nx.from_numpy_matrix(mat < m)
        else:
            G = nx.empty_graph(12)
            G.add_edges_from(edge_list)
            for (i, j) in ((i, j) for i in range(12) for j in range(12)):
                if ((i, j) not in edge_list) and ((j, i) not in edge_list):
                    mat[i][j] = m

        ax[0].imshow(mat, cmap='hot')
        ax[1].axis('off')
        
        nx.draw_circular(G, with_labels=True, ax=ax[1], edge_color=edge_color, width=edge_width)
    else:
        plt.imshow(mat, cmap='hot')
        plt.colorbar()
        plt.gcf().set_size_inches(6, 5)
    plt.savefig("../../presentation/figures/interaction_matrix_{}.pdf".format(suffix), transparent=True, dpi=200)
    return ax

if len(sys.argv) < 2:
    sys.argv.append("../../data/PNIC_shuffled.bin")

parser = argparse.ArgumentParser()
parser.add_argument('path', type=Path)
parser.add_argument('--nonbinding-cutoff', default=-7, type=float)
parser.add_argument('--binding-cutoff', default=-7.7, type=float)
args = parser.parse_args()

mat = load(str(args.path)).copy()
m = mat.max()

plot_graph(mat, m, "full", False)

mat[mat > args.nonbinding_cutoff] = m
ax = plot_graph(mat, m, "edges")

mat[mat > args.binding_cutoff] = m
plot_graph(mat, m, "strong_edges", ax=ax, edge_color='b', edge_width=2)

edge_list = [(2, 6), (3, 9), (4, 10), (7, 11)]
plot_graph(mat, m, "independent_edges", ax=ax, edge_color='g', edge_width=4, edge_list = edge_list)