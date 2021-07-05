import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager

from common import calculate_entropy, calculate_locational_variance, get_pdf

INPUT_ENTROPY_BEFORE, INPUT_ENTROPY_AFTER = snakemake.input.entropy
INPUT_VARIANCE_BEFORE, INPUT_VARIANCE_AFTER = snakemake.input.variance
INPUT_FONT_FILE = snakemake.input.font_file

print(snakemake.input.entropy)
print(snakemake.input.variance)

OUTPUT_OVERELAP_FIG = snakemake.output.overlap_ent_var_fig

max_entropy = snakemake.params.max_entropy
max_variance = snakemake.params.max_variance
n_entropy_bin = snakemake.params.n_entropy_bin
n_variance_bin = snakemake.params.n_variance_bin

entropy_before = np.load(INPUT_ENTROPY_BEFORE)
entropy_after = np.load(INPUT_ENTROPY_AFTER)
variance_before = np.load(INPUT_VARIANCE_BEFORE)
variance_after = np.load(INPUT_VARIANCE_AFTER)


# set bins
ent_bins = np.linspace(0, max_entropy, n_entropy_bin + 1)
var_bins = np.linspace(0, max_variance, n_variance_bin + 1)

# calculate pdf
ent_xs = [(ent_bins[i] + ent_bins[i + 1]) / 2 for i in range(len(ent_bins) - 1)]
var_xs = [(var_bins[i] + var_bins[i + 1]) / 2 for i in range(len(var_bins) - 1)]


ent_data_pdf_before = get_pdf(entropy_before, ent_bins)
ent_data_pdf_after = get_pdf(entropy_after, ent_bins)
var_data_pdf_before = get_pdf(variance_before, var_bins)
var_data_pdf_after = get_pdf(variance_after, var_bins)


# draw figure
plt.rcParams["figure.figsize"] = (12.0, 6.0)

prop = font_manager.FontProperties(fname=INPUT_FONT_FILE, size=22)
tiny_prop = font_manager.FontProperties(fname=INPUT_FONT_FILE, size=15)
fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.plot(
    ent_xs,
    ent_data_pdf_before,
    "-o",
    linewidth=1.5,
    markerfacecolor="white",
    markersize=0,
    color="tab:orange",
    alpha=0.7,
    label="Pre-COVID-19",
)
ax1.plot(
    ent_xs,
    ent_data_pdf_after,
    "-o",
    linewidth=1.5,
    markerfacecolor="white",
    markersize=0,
    color="tab:blue",
    alpha=0.7,
    label="Post-COVID-19",
)
ax1.set_xlabel("Hotspot entropy", fontproperties=prop)
ax1.set_ylabel("PDF", fontproperties=prop)

ax2.plot(
    var_xs,
    var_data_pdf_before,
    "-o",
    linewidth=1.5,
    markerfacecolor="white",
    markersize=0,
    color="tab:orange",
    alpha=0.7,
    label="Pre-COVID-19",
)
ax2.plot(
    var_xs,
    var_data_pdf_after,
    "-o",
    linewidth=1.5,
    markerfacecolor="white",
    markersize=0,
    color="tab:blue",
    alpha=0.7,
    label="Post-COVID-19",
)
ax2.set_xlabel("Radius of recreational activity", fontproperties=prop)

for ax in [ax1,ax2]:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

plt.legend(prop=prop, frameon=False)

for label in ax1.get_xticklabels():
    label.set_fontproperties(tiny_prop)
for label in ax1.get_yticklabels():
    label.set_fontproperties(tiny_prop)
for label in ax2.get_xticklabels():
    label.set_fontproperties(tiny_prop)
for label in ax2.get_yticklabels():
    label.set_fontproperties(tiny_prop)

plt.savefig(OUTPUT_OVERELAP_FIG, bbox_inches="tight")
