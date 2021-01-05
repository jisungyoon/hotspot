from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from matplotlib import cm, font_manager
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

INPUT_SEQUENCE = snakemake.input.sequence
INPUT_HOTSPOT_LEVEL_BY_GRID = snakemake.input.hotspot_level_by_grid
INPUT_VENDOR_TO_GRID = snakemake.input.vendor_to_grid
INPUT_FONT_FILE = snakemake.input.font_file

OUTPUT_HOTSPOT_MATRIX = snakemake.output.hotspot_matrix
OUTPUT_NULL_HOTSPOT_MATRIX = snakemake.output.null_hotspot_matrix

OUTPUT_HOTSPOT_MATRIX_FIG= snakemake.output.hotspot_matrix_fig
OUTPUT_NULL_HOTSPOT_MATRIX_FIG = snakemake.output.null_hotspot_matrix_fig

OUTPUT_PSI_RESULT = snakemake.output.psi_result

hotspot_level = snakemake.params.hotspot_level

# import data
sequences = np.load(INPUT_SEQUENCE, allow_pickle=True)
hotspot_level_by_grid = np.load(INPUT_HOTSPOT_LEVEL_BY_GRID)
vendor_to_grid = pd.read_pickle(INPUT_VENDOR_TO_GRID)

# construct hot spot transition matrix
hotspot_matrix = np.zeros((hotspot_level, hotspot_level))
for row in sequences:
    prev_level = int(hotspot_level_by_grid[vendor_to_grid[row[0]]])
    for x in row[1:]:
        cur_level = int(hotspot_level_by_grid[vendor_to_grid[x]])
        hotspot_matrix[prev_level - 1][cur_level - 1] += 1
        prev_level = cur_level
normed_hotspot_matrix = hotspot_matrix/np.sum(hotspot_matrix)
np.save(OUTPUT_HOTSPOT_MATRIX , normed_hotspot_matrix)

# draw hot spot trasition matrix
prop = font_manager.FontProperties(fname=INPUT_FONT_FILE, size=22)
small_prop = font_manager.FontProperties(fname=INPUT_FONT_FILE, size=20)
tiny_prop = font_manager.FontProperties(fname=INPUT_FONT_FILE, size=18)

f = plt.figure(figsize=(6.2,5.6))
ax = f.add_axes([0.17, 0.02, 0.72, 0.79])
axcolor = f.add_axes([0.93, 0.02, 0.03, 0.79])

im = ax.matshow(normed_hotspot_matrix, cmap=cm.Blues, norm=LogNorm())
cbar = f.colorbar(im, cax=axcolor)

ax.set_xticks(np.arange(0,10))
ax.set_xticklabels(np.arange(1,11))
ax.set_yticks(np.arange(0,10))
ax.set_yticklabels(np.arange(1,11))
ax.set_ylim((9.5, -0.5))

ax.set_ylabel('Hotspot level',fontproperties=prop)
ax.set_title('Hotspot level',fontproperties=prop, y=1.1)
for label in ax.get_xticklabels():
    label.set_fontproperties(small_prop)
for label in ax.get_yticklabels():
    label.set_fontproperties(small_prop)

cbar.ax.set_ylabel("Normalized flow", fontproperties=prop, labelpad =10)
for label in cbar.ax.get_yticklabels():
    label.set_fontproperties(tiny_prop)
plt.savefig(OUTPUT_HOTSPOT_MATRIX_FIG, bbox_inches='tight')

# calculate null hot spot matrix
null_hotspot_matrix = np.zeros((hotspot_level, hotspot_level))
total_sum = np.sum(hotspot_matrix)
for i, row in enumerate(null_hotspot_matrix):
    for j, x in enumerate(row):
        null_hotspot_matrix[i][j] = sum(hotspot_matrix[i, :]) * sum(hotspot_matrix[j, :]) / total_sum
normed_null_hotspot_matrix = null_hotspot_matrix/np.sum(null_hotspot_matrix)
np.save(OUTPUT_NULL_HOTSPOT_MATRIX , normed_null_hotspot_matrix)

# draw null hot spot trasition matrix
f = plt.figure(figsize=(6.2,5.6))
ax = f.add_axes([0.17, 0.02, 0.72, 0.79])
axcolor = f.add_axes([0.93, 0.02, 0.03, 0.79])

im = ax.matshow(normed_null_hotspot_matrix, cmap=cm.Blues, norm=LogNorm())
cbar = f.colorbar(im, cax=axcolor)

ax.set_xticks(np.arange(0,10))
ax.set_xticklabels(np.arange(1,11))
ax.set_yticks(np.arange(0,10))
ax.set_yticklabels(np.arange(1,11))
ax.set_ylim((9.5, -0.5))

ax.set_ylabel('Hotspot level',fontproperties=prop)
ax.set_title('Hotspot level',fontproperties=prop, y=1.1)
for label in ax.get_xticklabels():
    label.set_fontproperties(small_prop)
for label in ax.get_yticklabels():
    label.set_fontproperties(small_prop)

cbar.ax.set_ylabel("Normalized flow", fontproperties=prop, labelpad =10)
for label in cbar.ax.get_yticklabels():
    label.set_fontproperties(tiny_prop)
plt.savefig(OUTPUT_NULL_HOTSPOT_MATRIX_FIG, bbox_inches='tight')

def get_tridiagonal(A):
    ks = [-1, 0, 1]
    total_sum = 0
    for k in ks:
        total_sum += np.sum(np.diag(A, k=k))
    
    return total_sum

# caluclate psi and report
psi = get_tridiagonal(normed_hotspot_matrix)
psi_null = get_tridiagonal(normed_null_hotspot_matrix)
df = pd.DataFrame({
    'kind': ["data", "null_model"],
    'psi': [psi, psi_null]
})
df.to_csv(OUTPUT_PSI_RESULT, index=None)


