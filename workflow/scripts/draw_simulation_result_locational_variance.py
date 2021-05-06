import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from common import (calculate_entropy, calculate_locational_variance,
                    generate_sequence, get_pdf)

INPUT_BEFORE_FILES = snakemake.input.before_files
INPUT_AFTER_FILES = snakemake.input.after_files
FONT_PATH = snakemake.input.font_file

OUTPUT_SIMULATION_FIG_BY_GAMMA = snakemake.output.simulation_fig_by_gamma

gamms = snakemake.params.gammas


def get_results(files):
    ent_result_dict = {}
    var_result_dict = {}
    for file in files:
        temp_name = file.split('/')[-1].split('_')
        gamma = float(temp_name[1].split('.c')[0])

        temp = pd.read_csv(file)
        ent_result_dict[gamma] = (temp['ent_jsd'].mean(), temp['ent_jsd'].std())
        var_result_dict[gamma] = (temp['var_jsd'].mean(), temp['var_jsd'].std())
    
    return ent_result_dict, var_result_dict
    
    return ent_result_dict, var_result_dict

def find_best_results(result_dict):
    keys = list(result_dict.keys())
    val_list = np.array([result_dict[key][0] for key in keys])
    best_p, best_k = keys[np.argmin(val_list)]
    
    return best_p, best_k

before_ent_result_dict, before_var_result_dict = get_results(INPUT_BEFORE_FILES)
after_ent_result_dict, after_var_result_dict = get_results(INPUT_AFTER_FILES )

gammas = [np.round(x,2) for x in np.linspace(0, 5 , 201)]
before_vars = [ before_var_result_dict[gamma][0] for gamma in gammas]
after_vars = [ after_var_result_dict[gamma][0] for gamma in gammas]


plt.rcParams['figure.figsize'] = (8.0, 6.0)
font_path = FONT_PATH
prop = font_manager.FontProperties(fname=font_path, size=24)
small_prop = font_manager.FontProperties(fname=font_path, size=18)
tiny_prop = font_manager.FontProperties(fname=font_path, size=16)

fig, ax1 = plt.subplots(1, 1, sharey=True)

ax1.plot(gammas, before_vars, '--', markerfacecolor='white', color='darkorange', label='Pre-COVID-19')
ax1.plot(gammas, after_vars, '--', markerfacecolor='white', color='cornflowerblue', label='Post-COVID-19')
ax1.set_ylabel(r'Locational Variance JSD', fontproperties=prop)
fig.text(0.52, 0.02, r'$\gamma$', ha='center', fontproperties=prop)

plt.legend(fontsize=18, frameon=False)

for ax in [ax1]:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    for label in ax.get_xticklabels():
        label.set_fontproperties(tiny_prop)
            
    for label in ax.get_yticklabels():
        label.set_fontproperties(tiny_prop)


plt.savefig(OUTPUT_SIMULATION_FIG_BY_GAMMA, bbox_inches='tight')