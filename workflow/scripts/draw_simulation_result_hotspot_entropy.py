import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (calculate_entropy, calculate_locational_variance,
                    generate_sequence, get_pdf)

INPUT_BEFORE_FILES = snakemake.input.before_files
INPUT_AFTER_FILES = snakemake.input.after_files
FONT_PATH = snakemake.input.font_file

OUTPUT_SIMULATION_FIG_BY_PK = snakemake.output.simulation_fig_by_pk

ps = snakemake.params.ps
ks = snakemake.params.ks


def get_results(files):
    ent_result_dict = {}
    var_result_dict = {}
    for file in files:
        temp_name = file.split("/")[-1].split("_")
        p = float(temp_name[1])
        k = float(temp_name[2].split(".c")[0])

        temp = pd.read_csv(file)
        ent_result_dict[(p, k)] = (temp["ent_jsd"].mean(), temp["ent_jsd"].std())
        var_result_dict[(p, k)] = (temp["var_jsd"].mean(), temp["var_jsd"].std())

    return ent_result_dict, var_result_dict


def find_best_results(result_dict):
    keys = list(result_dict.keys())
    val_list = np.array([result_dict[key][0] for key in keys])
    best_p, best_k = keys[np.argmin(val_list)]

    return best_p, best_k


before_ent_result_dict, before_var_result_dict = get_results(INPUT_BEFORE_FILES)
after_ent_result_dict, after_var_result_dict = get_results(INPUT_AFTER_FILES)

before_best_p, before_best_k = find_best_results(before_ent_result_dict)
after_best_p, after_best_k = find_best_results(after_ent_result_dict)

ent_jsds_before = np.zeros((len(ps), len(ks)))
for i, p in enumerate(ps):
    for j, k in enumerate(ks):
        ent_jsds_before[i][j] = before_ent_result_dict[(p, k)][0]

ent_jsds_after = np.zeros((len(ps), len(ks)))
for i, p in enumerate(ps):
    for j, k in enumerate(ks):
        ent_jsds_after[i][j] = after_ent_result_dict[(p, k)][0]

plt.rcParams["figure.figsize"] = (13.0, 6.0)

prop = font_manager.FontProperties(fname=FONT_PATH, size=25)
small_prop = font_manager.FontProperties(fname=FONT_PATH, size=18)
tiny_prop = font_manager.FontProperties(fname=FONT_PATH, size=18)

fig, (ax1, ax2) = plt.subplots(1, 2)

im = ax1.pcolormesh(ent_jsds_before, vmin=0.053, vmax=0.531)
im_ = ax2.pcolormesh(ent_jsds_after, vmin=0.053, vmax=0.531)
for ax in ax1, ax2:
    ax.set_yticks([0, 11, 21, 31, 41, 51])
    ax.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontproperties=tiny_prop)
    ax.set_xticks([0, 20, 40, 60, 80, 100, 121])
    ax.set_xticklabels(
        [ks[0], ks[20], ks[40], ks[60], ks[80], ks[100], ks[120]],
        fontproperties=tiny_prop,
    )
ax1.set_ylabel(r"$p$", fontproperties=prop)
fig.text(0.465, 0.02, r"$k$", ha="center", fontproperties=prop)

min_p, min_k = np.unravel_index(ent_jsds_before.argmin(), ent_jsds_before.shape)
ax1.plot(min_k, min_p, "r+", markersize=20)
ax1.text(62, 7.5, r"$p^*={0:.3f}$".format(before_best_p), fontproperties=prop)
ax1.text(62, 2, r"$k^*={0:.3f}$".format(before_best_k), fontproperties=prop)

min_p, min_k = np.unravel_index(ent_jsds_after.argmin(), ent_jsds_after.shape)
ax2.plot(min_k, min_p, "r+", markersize=20)
ax2.text(62, 7.5, r"$p^*={0:.3f}$".format(after_best_p), fontproperties=prop)
ax2.text(62, 2, r"$k^*={0:.3f}$".format(after_best_k), fontproperties=prop)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.83, 0.125, 0.015, 0.755])
cbar = fig.colorbar(im, cax=cbar_ax)


for label in cbar.ax.get_yticklabels():
    label.set_fontproperties(tiny_prop)
cbar.ax.set_ylabel("Hotspot Entropy JSD", fontproperties=prop, labelpad=13)
plt.savefig(OUTPUT_SIMULATION_FIG_BY_PK, bbox_inches="tight")
