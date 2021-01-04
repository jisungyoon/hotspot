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

###############################################################################
# DERIVED_DATA
###############################################################################
DERIVED_DATA_BY_PERIOD = j(DERIVED_DIR, "{period}")
SEQUENCE_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "sequence.npy")
SEQUENCE_LENGTH_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "sequence_length.npy")


###############################################################################
# FIGS
###############################################################################




rule all:
    input:
        expand(SEQUENCE_BY_PERIOD, period=PERIODS),
        expand(SEQUENCE_LENGTH_BY_PERIOD, period=PERIODS)


rule generate_sequence:
    input: RAW_SEQUENCE_DATA
    output: sequence=SEQUENCE_BY_PERIOD, sequence_length=SEQUENCE_LENGTH_BY_PERIOD
    script: "workflow/scripts/generate_sequence.py"