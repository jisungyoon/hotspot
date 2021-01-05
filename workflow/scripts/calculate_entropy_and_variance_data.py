import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

from common import calculate_entropy, calculate_locational_variance, get_pdf

INPUT_SEQUENCE = snakemake.input.sequence
INPUT_VENDOR_TO_GRID= snakemake.input.vendor_to_grid
INPUT_HOTSPOT_LEVEL_BY_GRID = snakemake.input.hotspot_level_by_grid
INPUT_FONT_FILE = snakemake.input.font_file

OUTPUT_ENTROPY = snakemake.output.entropy
OUTPUT_VARAINCE = snakemake.output.variance
OUTPUT_ENT_VAR_FIG = snakemake.output.ent_var_fig

hotspot_level = snakemake.params.hotspot_level

max_entropy = snakemake.params.max_entropy
max_variance = snakemake.params.max_variance
n_entropy_bin = snakemake.params.n_entropy_bin
n_variance_bin = snakemake.params.n_variance_bin

sequences = np.load(INPUT_SEQUENCE, allow_pickle=True)
vendor_to_grid = pd.read_pickle(INPUT_VENDOR_TO_GRID)
hotspot_level_by_grid = np.load(INPUT_HOTSPOT_LEVEL_BY_GRID)
grid_to_hot_spot_level = {(i,j): val for i, row in enumerate(hotspot_level_by_grid) for j, val in enumerate(row) if val != 11}

sequences = [list(map(vendor_to_grid.get, row)) for row in sequences]

# calculate entropy and variance
entropies = calculate_entropy(sequences, grid_to_hot_spot_level, hotspot_level=hotspot_level)
variances = np.mean([calculate_locational_variance(sequences) for i in range(10)], axis=0)

# save result
np.save(OUTPUT_ENTROPY, entropies)
np.save(OUTPUT_VARAINCE, variances)

# set bins
ent_bins = np.linspace(0, max_entropy, n_entropy_bin + 1)
var_bins = np.linspace(0, max_variance, n_variance_bin + 1)

# calculate pdf
ent_xs = [(ent_bins[i] + ent_bins[i + 1])/2 for i in range(len(ent_bins) - 1)]
var_xs = [(var_bins[i] + var_bins[i + 1])/2 for i in range(len(var_bins) - 1)]
var_data_pdf = get_pdf(variances, var_bins)
ent_data_pdf = get_pdf(entropies, ent_bins)

# draw figure
plt.rcParams['figure.figsize'] = (12.0, 6.0)

prop = font_manager.FontProperties(fname=INPUT_FONT_FILE , size=22)
tiny_prop = font_manager.FontProperties(fname=INPUT_FONT_FILE, size=15)
fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.plot(ent_xs, ent_data_pdf, '-o', markerfacecolor='white', color='darkorange', label='data')
ax1.set_xlabel('Entropy', fontproperties=prop)
ax1.set_ylabel('PDF', fontproperties=prop)

ax2.plot(var_xs, var_data_pdf, '-o', markerfacecolor='white', color='darkorange', label='data')
ax2.set_xlabel('Locational Variance', fontproperties=prop)

plt.legend(prop=prop, frameon=False)

for label in ax1.get_xticklabels():
    label.set_fontproperties(tiny_prop)
for label in ax1.get_yticklabels():
    label.set_fontproperties(tiny_prop)
for label in ax2.get_xticklabels():
    label.set_fontproperties(tiny_prop)
for label in ax2.get_yticklabels():
    label.set_fontproperties(tiny_prop)
    
plt.savefig(OUTPUT_ENT_VAR_FIG, bbox_inches='tight')