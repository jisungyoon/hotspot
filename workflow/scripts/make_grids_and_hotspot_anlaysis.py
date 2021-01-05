from collections import Counter, defaultdict

import numpy as np
import pandas as pd

INPUT_RAW = snakemake.input.raw
INPUT_META = snakemake.input.meta

OUTPUT_RESERVATION_BY_GRID = snakemake.output.reservation_by_grid
OUTPUT_HOT_SPOT_LEVEL_BY_GRID = snakemake.output.hot_spot_level_by_grid

OUTPUT_VENDOR_TO_GRID = snakemake.output.vendor_to_grid
OUTPUT_VENDOR_TO_HOT_SPOT_LEVEL = snakemake.output.vendor_to_hotspot_level


n_lat_bin = snakemake.params.n_lat_bin
n_lon_bin = snakemake.params.n_lon_bin

max_lat = snakemake.params.max_lat
min_lat = snakemake.params.min_lat
max_lon = snakemake.params.max_lon
min_lon = snakemake.params.min_lon

hot_spot_level = snakemake.params.hot_spot_level


# Import data
data = pd.read_csv(INPUT_RAW, sep="|", parse_dates=["dt_reservation", "dt_entrance"])
data = data[data.order_cnt > 0]  # using only valid reservation
data = data[data.uno2 != 0.0]  # remove Non-Members

# Import meta_data
meta_data = pd.read_csv(INPUT_META, sep="|",).set_index("ano")

# get locational info
vendors = set(data.yg_vendor_code)
lat_lon_dict = {}
for id_, lon, lat in zip(meta_data.index, meta_data.alng, meta_data.alat):
    lat_lon_dict[id_] = (lat, lon)
lats = [lat_lon_dict[vendor][0] for vendor in vendors]
lons = [lat_lon_dict[vendor][1] for vendor in vendors]


# make grids
lat_bins = np.linspace(min_lat, max_lat, n_lat_bin + 1)
lon_bins = np.linspace(min_lon, max_lon, n_lon_bin + 1)
lat_digitized = np.digitize(lats, lat_bins) - 1
lon_digitized = np.digitize(lons, lon_bins) - 1

# count reservation by grid
reservation_by_grid = np.zeros((n_lon_bin, n_lat_bin))
frequency_counter = data.yg_vendor_code.value_counts()
for i, vendor in enumerate(vendors):
    lat_bin = lat_digitized[i]
    lon_bin = lon_digitized[i]
    reservation_by_grid[lon_bin][lat_bin] += frequency_counter[vendor]
np.save(OUTPUT_RESERVATION_BY_GRID, reservation_by_grid)

# assign hot spot level
non_zero_values = reservation_by_grid[reservation_by_grid != 0]
hot_spot_level_by_grid = np.zeros((n_lon_bin, n_lat_bin))
percentiles = np.linspace(0, 100, hot_spot_level + 1)
for i, per in enumerate(percentiles):
    th = np.percentile(non_zero_values, per)
    if i == 0:
        th -= 1
    hot_spot_level_by_grid[reservation_by_grid > th] = i + 1
hot_spot_level_by_grid = 11 - hot_spot_level_by_grid
np.save(OUTPUT_HOT_SPOT_LEVEL_BY_GRID, hot_spot_level_by_grid)


# vendor matching data files
vendor_to_grid = {
    vendor: (lon_digitized[i], lat_digitized[i]) for i, vendor in enumerate(vendors)
}
vendor_to_hot_spot_level = {
    vendor: hot_spot_level_by_grid[vendor_to_grid[vendor]] for vendor in vendors
}
pd.to_pickle(vendor_to_grid, OUTPUT_VENDOR_TO_GRID)
pd.to_pickle(vendor_to_hot_spot_level, OUTPUT_VENDOR_TO_HOT_SPOT_LEVEL)
