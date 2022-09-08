from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon

from common import (
    calculate_entropy,
    calculate_locational_variance,
    calculate_p,
    generate_sequence,
    get_pdf,
)

INPUT_HOME_PDF = snakemake.input.home_pdf
INPUT_SEQUENCE_LENGTH = snakemake.input.sequence_length
INPUT_HOTSPOT_LEVEL_BY_GRID = snakemake.input.hotspot_level_by_grid
INPUT_DISTANCE = snakemake.input.dist
INPUT_ENTROPY_DATA = snakemake.input.entropy_data
INPUT_P_DATA = snakemake.input.p_data
INPUT_VARIANCE_DATA = snakemake.input.variance_data

OUTPUT_JSD_RESULT = snakemake.output.jsd
OUTPUT_ENTROPY = snakemake.output.entropy
OUTPUT_VARAINCE = snakemake.output.variance

hotspot_level = snakemake.params.hotspot_level

max_entropy = snakemake.params.max_entropy
max_variance = snakemake.params.max_variance
n_entropy_bin = snakemake.params.n_entropy_bin
n_variance_bin = snakemake.params.n_variance_bin

p = float(snakemake.wildcards.p)
k = float(snakemake.wildcards.k)
gamma = 1

repetition = int(snakemake.params.repetition)


home_pdf = pd.read_pickle(INPUT_HOME_PDF)
sequence_length = np.load(INPUT_SEQUENCE_LENGTH, allow_pickle=True)

grid_to_hotspot_level = pd.read_pickle(INPUT_HOTSPOT_LEVEL_BY_GRID)
hotspot_level_to_grid = defaultdict(list)
for k_, v_ in grid_to_hotspot_level.items():
    hotspot_level_to_grid[v_].append(k_)

homes = list(home_pdf.keys())
homes_p = [home_pdf[home] for home in homes]

d = pd.read_pickle(INPUT_DISTANCE)

entropies_array = []
variances_array = []
ps_array = []

# generate sequences
for i in range(repetition):
    generated_homes_idx = np.random.choice(
        np.arange(len(homes)), size=len(sequence_length), replace=True, p=homes_p
    )
    generated_homes = [homes[idx] for idx in generated_homes_idx]
    generated_sequences = generate_sequence(
        p,
        gamma,
        k,
        generated_homes,
        sequence_length,
        hotspot_level_to_grid,
        d,
    )

    entropies = calculate_entropy(
        generated_sequences, grid_to_hotspot_level, hotspot_level=hotspot_level
    )
    variances = np.mean(
        [calculate_locational_variance(generated_sequences, d) for i in range(10)],
        axis=0,
    )
    ps = calculate_p(generated_sequences, grid_to_hotspot_level, hotspot_level=10)

    entropies_array.append(entropies)
    variances_array.append(variances)
    ps_array.append(ps)

# calculate JSD
entropy_data = np.load(INPUT_ENTROPY_DATA)
variance_data = np.load(INPUT_VARIANCE_DATA)
ps_data = np.load(INPUT_P_DATA)

# set bins
ent_bins = np.linspace(0, max_entropy, n_entropy_bin + 1)
var_bins = np.linspace(0, max_variance, n_variance_bin + 1)

# calculate data pdf
ent_data_pdf = get_pdf(entropy_data, ent_bins)
var_data_pdf = get_pdf(variance_data, var_bins)

ent_jsds = []
var_jsds = []
ps_jsds = []

for entropies, variances, ps in zip(entropies_array, variances_array, ps_array):
    ent_sim_pdf = get_pdf(entropies, ent_bins)
    var_sim_pdf = get_pdf(variances, var_bins)

    ent_jsd = jensenshannon(ent_data_pdf, ent_sim_pdf)
    var_jsd = jensenshannon(var_data_pdf, var_sim_pdf)
    ps_jsd = jensenshannon(ps_data, ps)

    ent_jsds.append(ent_jsd)
    var_jsds.append(var_jsd)
    ps_jsds.append(ps_jsd)


df = pd.DataFrame(
    {
        "ent_jsd": ent_jsds,
        "var_jsd": var_jsds,
        "p_jsd": ps_jsd,
        "repetition": list(range(repetition)),
    }
)
df.to_csv(OUTPUT_JSD_RESULT, index=None)
np.save(OUTPUT_ENTROPY, entropies_array)
np.save(OUTPUT_VARAINCE, variances_array)
