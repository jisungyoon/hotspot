from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import s2geometry as s2

INPUT_RAW = snakemake.input.raw
INPUT_META = snakemake.input.meta

OUTPUT_RESERVATION_BY_GRID = snakemake.output.reservation_by_grid
OUTPUT_HOTSPOT_LEVEL_BY_GRID = snakemake.output.hotspot_level_by_grid

OUTPUT_VENDOR_TO_GRID = snakemake.output.vendor_to_grid
OUTPUT_VENDOR_TO_HOTSPOT_LEVEL = snakemake.output.vendor_to_hotspot_level

hotspot_level = snakemake.params.hotspot_level


# Import data
data = pd.read_csv(INPUT_RAW, sep="|", parse_dates=["dt_reservation", "dt_entrance"])
data = data[data.order_cnt > 0]  # using only valid reservation
data = data[data.uno2 != 0.0]  # remove Non-Members

# Import meta_data
meta_data = pd.read_csv(INPUT_META, sep="|",).set_index("ano")

# get locational info
vendors = set(data.yg_vendor_code)
meta_data = meta_data[meta_data.index.isin(vendors)]
vendor_to_grid = {}
for id_, lon, lat in zip(meta_data.index, meta_data.alng, meta_data.alat):
    vendor_to_grid[id_] = s2.S2CellId(s2.S2LatLng.FromDegrees(lat, lon)).parent(14).id()
pd.to_pickle(vendor_to_grid, OUTPUT_VENDOR_TO_GRID)

# count reservation by grid
reservation_counter = data.yg_vendor_code.value_counts().to_dict()
reservation_counter_by_grid = defaultdict(int)
for k, v in reservation_counter.items():
    reservation_counter_by_grid[vendor_to_grid[k]] += v
pd.to_pickle(reservation_counter_by_grid, OUTPUT_RESERVATION_BY_GRID)

# assign hot spot level
hot_spot_level = 10
ks = np.array(list(reservation_counter_by_grid.keys()))
vals = np.array(list(reservation_counter_by_grid.values()))

hotspot_level_by_grid = {}
percentiles = np.linspace(0, 100, hot_spot_level + 1)
for i, per in enumerate(percentiles):
    th = np.percentile(vals, per)
    if i == 0:
        th -= 1
    for k in ks[vals > th]:
        hotspot_level_by_grid[k] = 10 - i
pd.to_pickle(hotspot_level_by_grid, OUTPUT_HOTSPOT_LEVEL_BY_GRID)


# vendor matching data files
vendor_to_hotspot_level = {
    vendor: hotspot_level_by_grid[vendor_to_grid[vendor]] for vendor in vendors
}

pd.to_pickle(vendor_to_hotspot_level, OUTPUT_VENDOR_TO_HOTSPOT_LEVEL)
