import numpy as np
import pandas as pd

from common import calculate_p

INPUT_SEQUENCE = snakemake.input.sequence
INPUT_MASKING = snakemake.input.masking
INPUT_VENDOR_TO_GRID = snakemake.input.vendor_to_grid
INPUT_GRID_TO_HOTSPOT_LEVEL = snakemake.input.grid_to_hotspot_level

OUTPUT_P = snakemake.output.p

hotspot_level = snakemake.params.hotspot_level

sequences = np.load(INPUT_SEQUENCE, allow_pickle=True)
masking = np.load(INPUT_MASKING)
sequences = sequences[masking]
vendor_to_grid = pd.read_pickle(INPUT_VENDOR_TO_GRID)
grid_to_hotspot_level = pd.read_pickle(INPUT_GRID_TO_HOTSPOT_LEVEL)

sequences = [list(map(vendor_to_grid.get, row)) for row in sequences]

# calculate entropy and variance
ps = calculate_p(sequences, grid_to_hotspot_level, hotspot_level=hotspot_level)
np.save(OUTPUT_P, ps)
