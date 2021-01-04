from os.path import join as j

configfile: "workflow/config.yaml"

DATA_DIR = config["data_dir"]
RAW_DIR = j(DATA_DIR, "raw")
DERIVED_DIR = j(DATA_DIR, "derived")
SHAPE_FILE_DIR = j(DATA_DIR, "shape_file", "seoul")

FIGURE_DIR = config["figs_dir"]
ASSETS_DIR = config["assets_dir"]

###############################################################################
# RAW_DATA
###############################################################################
PERIODS = ["before", "after"]
RAW_SEQUENCE_DATA = j(RAW_DIR, "{period}_data.csv")
META_DATA = j(RAW_DIR, "meta_data.csv")

###############################################################################
# DERIVED_DATA
###############################################################################
DERIVED_DATA_BY_PERIOD = j(DERIVED_DIR, "{period}")
SEQUENCE_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "sequence.npy")
SEQUENCE_LENGTH_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "sequence_length.npy")

GRID_RESERVATION_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "reservation_by_grid.npy")
GRID_HOT_SPOT_LEVEL_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "hot_spot_level_by_grid.npy")

VENDOR_TO_GRID_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "vendor_to_grid.pkl")
VENDOR_TO_HOT_SPOT_LEVEL_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "vendor_to_hot_spot_level.pkl")

###############################################################################
# FIGS
###############################################################################
FIGS_BY_PERIOD = j(FIGURE_DIR, "{period}")



###############################################################################
# PARAMS
###############################################################################
N_LAT_BIN = 30
N_LON_BIN = 30

MAX_LAT = 37.69
MIN_LAT = 37.45
MAX_LON = 127.16
MIN_LON = 126.80

HOT_SPOT_LEVEL = 10





rule all:
    input:
        expand(SEQUENCE_BY_PERIOD, period=PERIODS),
        expand(SEQUENCE_LENGTH_BY_PERIOD, period=PERIODS),
        expand(GRID_RESERVATION_BY_PERIOD, period=PERIODS),
        expand(GRID_HOT_SPOT_LEVEL_BY_PERIOD, period=PERIODS),
        expand(VENDOR_TO_GRID_BY_PERIOD, period=PERIODS),
        expand(VENDOR_TO_HOT_SPOT_LEVEL_BY_PERIOD, period=PERIODS)


rule generate_sequence:
    input: RAW_SEQUENCE_DATA
    output: sequence=SEQUENCE_BY_PERIOD, sequence_length=SEQUENCE_LENGTH_BY_PERIOD
    script: "workflow/scripts/generate_sequence.py"
    
    
rule hotspot_analysis:
    input: raw=RAW_SEQUENCE_DATA, meta=META_DATA
    output: reservation_by_grid=GRID_RESERVATION_BY_PERIOD, hot_spot_level_by_grid=GRID_HOT_SPOT_LEVEL_BY_PERIOD, vendor_to_grid=VENDOR_TO_GRID_BY_PERIOD, vendor_to_hotspot_level=VENDOR_TO_HOT_SPOT_LEVEL_BY_PERIOD
    params: n_lat_bin=N_LAT_BIN, n_lon_bin = N_LON_BIN, max_lat=MAX_LAT, min_lat=MIN_LAT, max_lon=MAX_LON, min_lon=MIN_LON, hot_spot_level=HOT_SPOT_LEVEL 
    script: "workflow/scripts/make_grids_and_hotspot_anlaysis.py"