from collections import Counter, defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm, figure, font_manager
from matplotlib.colors import LogNorm
from scipy.stats import pearsonr

INPUT_SEQUENCE = snakemake.input.sequence
INPUT_VENDOR_TO_GRID = snakemake.input.vendor_to_grid
INPUT_RESERVATION_BY_GRID = snakemake.input.reservation_by_grid
INPUT_FONT_FILE = snakemake.input.font_file

OUTPUT_HOME_PDF = snakemake.output.home_pdf
OUTPUT_HOME_VALIDATION_FIG = snakemake.output.home_validation_fig

sequences = np.load(INPUT_SEQUENCE, allow_pickle=True)
vendor_to_grid = pd.read_pickle(INPUT_VENDOR_TO_GRID)
reservation_by_grid = pd.read_pickle(INPUT_RESERVATION_BY_GRID)

home_list = []
for row in sequences:
    counter_dict = dict(Counter(list(map(vendor_to_grid.get, row))))
    max_val = max(counter_dict.values())
    max_keys = [k for k, v in counter_dict.items() if v == max_val]
    home = max_keys[np.random.choice(len(max_keys))]
    home_list.append(home)

# valdiate hypothesis
home_counter = dict(Counter(home_list))
xs = [v for k, v in home_counter.items()]
ys = [reservation_by_grid[k] for k, v in home_counter.items()]

prop = font_manager.FontProperties(fname=INPUT_FONT_FILE, size=22)
tiny_prop = font_manager.FontProperties(fname=INPUT_FONT_FILE, size=15)

plt.hexbin(
    xs, ys, xscale="log", yscale="log", gridsize=(40, 15), mincnt=1, cmap="Greys"
)
plt.ylabel(r"$N_r$", fontproperties=prop)
plt.xlabel(r"$N_h$", fontproperties=prop)
plt.xticks(fontproperties=tiny_prop)
plt.yticks(fontproperties=tiny_prop)
plt.xlim([2, 15000])
plt.ylim([10, 200000])

plt.text(
    3, 60000, "r = {}".format(np.round(pearsonr(xs, ys)[0], 3)), fontproperties=prop
)

cbar = plt.colorbar()
cbar.ax.set_ylabel("Count", fontproperties=prop)
for label in cbar.ax.get_yticklabels():
    label.set_fontproperties(tiny_prop)
plt.xlim([2, 15000])
plt.ylim([10, 200000])
plt.savefig(OUTPUT_HOME_VALIDATION_FIG, bbox_inches="tight")


# calculate pdf
homes_p_dict = {}
for k, v in reservation_by_grid.items():
    if v > 0:
        homes_p_dict[k] = v

total_sum = sum(homes_p_dict.values())
for k, v in homes_p_dict.items():
    homes_p_dict[k] = v / total_sum
pd.to_pickle(homes_p_dict, OUTPUT_HOME_PDF)
