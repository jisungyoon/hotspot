from os.path import join as j
import numpy as np

configfile: "workflow/config.yaml"

DATA_DIR = config["data_dir"]
SIMULATION_DIR = config["simulation_dir"]
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
SHAPE_FILE = j(SHAPE_FILE_DIR , "seoul.shp")

###############################################################################
# DERIVED_DATA
###############################################################################
DERIVED_DATA_BY_PERIOD = j(DERIVED_DIR, "{period}")
SEQUENCE_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "sequence.npy")
SEQUENCE_LENGTH_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "sequence_length.npy")

GRID_RESERVATION_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "reservation_by_grid.npy")
GRID_HOTSPOT_LEVEL_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "hotspot_level_by_grid.npy")

VENDOR_TO_GRID_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "vendor_to_grid.pkl")
VENDOR_TO_HOTSPOT_LEVEL_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "vendor_to_hotspot_level.pkl")

HOTSPOT_MATRIX_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "hotspot_matrix.npy")
NULL_HOTSPOT_MATRIX_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "null_hotspot_matrix.npy")
PSI_RESULT_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "psi.csv")

HOME_PDF_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "home_pdf.pkl")

ENTROPY_DATA_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "entropy_data.npy")
VARIANCE_DATA_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD , "varaince_data.npy")

###############################################################################
# SIMULATION_RESULT
###############################################################################
SIMULATION_RESULT_DIR = j(SIMULATION_DIR, "{period}")

SIMULATION_JSD_RESULT = j(SIMULATION_RESULT_DIR, "result")
SIMULATION_ENTROPY_RESULT= j(SIMULATION_RESULT_DIR, "entropy")
SIMULATION_VARIANCE_RESULT= j(SIMULATION_RESULT_DIR, "variance")

SIMULATION_JSD_RESULT_FILE = j(SIMULATION_JSD_RESULT, "result_{p}_{k}.csv")
SIMULATION_ENTROPY_RESULT_FILE = j(SIMULATION_ENTROPY_RESULT, "entropy_{p}_{k}.npy")
SIMULATION_VARIANCE_RESULT_FILE = j(SIMULATION_VARIANCE_RESULT, "variance_{p}_{k}.npy")


###############################################################################
# ASSETS
###############################################################################
NORMAL_FONT_PATH = j(ASSETS_DIR, "Helvetica.ttf")

###############################################################################
# FIGS
###############################################################################
FIGS_BY_PERIOD = j(FIGURE_DIR, "{period}")

RESERVATION_MAP_BY_PERIOD = j(FIGS_BY_PERIOD, "reservation_map.pdf")
HOTSPOT_MAP_BY_PERIOD = j(FIGS_BY_PERIOD, "hotspot_map.pdf")

HOTSPOT_MATRIX_FIG_BY_PERIOD = j(FIGS_BY_PERIOD, "hotspot_matrix.pdf")
NULL_HOTSPOT_MATRIX_FIG_BY_PERIOD = j(FIGS_BY_PERIOD, "null_hotspot_matrix.pdf")

HOME_VALIDATION_BY_PERIOD = j(FIGS_BY_PERIOD, "home_validation.pdf")

ENT_VAR_BY_PERIOD = j(FIGS_BY_PERIOD, "entropy_variance_data.pdf")
OVERLAP_ENT_VAR = j(FIGURE_DIR, "entropy_variance_data_overlap.pdf")


###############################################################################
# PARAMS
###############################################################################
N_LAT_BIN = 30
N_LON_BIN = 30

MAX_LAT = 37.69
MIN_LAT = 37.45
MAX_LON = 127.16
MIN_LON = 126.80

HOTSPOT_LEVEL = 10

MAX_ENTROPY = 2.2
MAX_VARIANCE = 24
 
N_ENTROPY_BIN = 30
N_VARIANCE_BIN = 30

###############################################################################
# PARAMS for simulation
###############################################################################
ps = [np.round(x,2) for x in np.linspace(0,1, 101)]
ks = [np.round(x,2) for x in np.linspace(0.05, 5 , 100)]
REPTITION = 10


rule all:
    input:
        expand(SEQUENCE_BY_PERIOD, period=PERIODS),
        expand(SEQUENCE_LENGTH_BY_PERIOD, period=PERIODS),
        expand(RESERVATION_MAP_BY_PERIOD, period=PERIODS),
        expand(HOTSPOT_MAP_BY_PERIOD, period=PERIODS),
        expand(HOTSPOT_MATRIX_FIG_BY_PERIOD, period=PERIODS),
        expand(NULL_HOTSPOT_MATRIX_FIG_BY_PERIOD, period=PERIODS),
        expand(PSI_RESULT_BY_PERIOD, period=PERIODS),
        expand(HOME_VALIDATION_BY_PERIOD, period=PERIODS),
        expand(ENT_VAR_BY_PERIOD, period=PERIODS),
        OVERLAP_ENT_VAR,
        expand(SIMULATION_JSD_RESULT_FILE, period=PERIODS, p=ps, k=ks),
        expand(SIMULATION_ENTROPY_RESULT_FILE, period=PERIODS, p=ps, k=ks),
        expand(SIMULATION_VARIANCE_RESULT_FILE, period=PERIODS, p=ps, k=ks)
      
rule generate_sequence:
    input: RAW_SEQUENCE_DATA
    output: sequence=SEQUENCE_BY_PERIOD, sequence_length=SEQUENCE_LENGTH_BY_PERIOD
    script: "workflow/scripts/generate_sequence.py"
    
    
rule hotspot_analysis:
    input: raw=RAW_SEQUENCE_DATA, meta=META_DATA
    output: reservation_by_grid=GRID_RESERVATION_BY_PERIOD, hotspot_level_by_grid=GRID_HOTSPOT_LEVEL_BY_PERIOD, vendor_to_grid=VENDOR_TO_GRID_BY_PERIOD, vendor_to_hotspot_level=VENDOR_TO_HOTSPOT_LEVEL_BY_PERIOD
    params: n_lat_bin=N_LAT_BIN, n_lon_bin = N_LON_BIN, max_lat=MAX_LAT, min_lat=MIN_LAT, max_lon=MAX_LON, min_lon=MIN_LON, hotspot_level=HOTSPOT_LEVEL 
    script: "workflow/scripts/make_grids_and_hotspot_anlaysis.py"
    
rule draw_reservation_map:
    input: reservation_by_grid=GRID_RESERVATION_BY_PERIOD, hotspot_level_by_grid=GRID_HOTSPOT_LEVEL_BY_PERIOD, shape_file=SHAPE_FILE, font_file=NORMAL_FONT_PATH
    output: reservation_map=RESERVATION_MAP_BY_PERIOD, hotspot_map=HOTSPOT_MAP_BY_PERIOD
    params: n_lat_bin=N_LAT_BIN, n_lon_bin = N_LON_BIN, max_lat=MAX_LAT, min_lat=MIN_LAT, max_lon=MAX_LON, min_lon=MIN_LON
    script: "workflow/scripts/draw_reservation_map.py"
    
    
rule calculate_transition_matrix_and_psi:
    input: sequence=SEQUENCE_BY_PERIOD, hotspot_level_by_grid=GRID_HOTSPOT_LEVEL_BY_PERIOD, vendor_to_grid=VENDOR_TO_GRID_BY_PERIOD, font_file=NORMAL_FONT_PATH
    output: hotspot_matrix=HOTSPOT_MATRIX_BY_PERIOD, null_hotspot_matrix=NULL_HOTSPOT_MATRIX_BY_PERIOD, hotspot_matrix_fig=HOTSPOT_MATRIX_FIG_BY_PERIOD, null_hotspot_matrix_fig=NULL_HOTSPOT_MATRIX_FIG_BY_PERIOD, psi_result=PSI_RESULT_BY_PERIOD
    params: hotspot_level=HOTSPOT_LEVEL 
    script: "workflow/scripts/calculate_transition_matrix_and_psi.py"
    
    
rule calculate_and_validate_home_distribution:
    input: sequence=SEQUENCE_BY_PERIOD, vendor_to_grid=VENDOR_TO_GRID_BY_PERIOD, reservation_by_grid=GRID_RESERVATION_BY_PERIOD, font_file=NORMAL_FONT_PATH
    output: home_pdf=HOME_PDF_BY_PERIOD, home_validation_fig=HOME_VALIDATION_BY_PERIOD 
    script: "workflow/scripts/calculate_and_validate_home_distribution.py"

rule calculate_entropy_variance_data:
    input: sequence=SEQUENCE_BY_PERIOD, vendor_to_grid=VENDOR_TO_GRID_BY_PERIOD, hotspot_level_by_grid=GRID_HOTSPOT_LEVEL_BY_PERIOD, font_file=NORMAL_FONT_PATH
    output: entropy=ENTROPY_DATA_BY_PERIOD, variance=VARIANCE_DATA_BY_PERIOD, ent_var_fig=ENT_VAR_BY_PERIOD
    params: hotspot_level=HOTSPOT_LEVEL, max_entropy=MAX_ENTROPY, max_variance=MAX_VARIANCE, n_entropy_bin=N_ENTROPY_BIN, n_variance_bin=N_VARIANCE_BIN
    script: "workflow/scripts/calculate_entropy_and_variance_data.py"
    
rule draw_overlap_ent_var_fig:
    input: entropy=expand(ENTROPY_DATA_BY_PERIOD, period=PERIODS), variance=expand(VARIANCE_DATA_BY_PERIOD, period=PERIODS), font_file=NORMAL_FONT_PATH
    output: overlap_ent_var_fig=OVERLAP_ENT_VAR 
    params: hotspot_level=HOTSPOT_LEVEL, max_entropy=MAX_ENTROPY, max_variance=MAX_VARIANCE, n_entropy_bin=N_ENTROPY_BIN, n_variance_bin=N_VARIANCE_BIN
    script: "workflow/scripts/draw_overlapped_ent_var.py"
    

rule simulation:
    input: home_pdf=HOME_PDF_BY_PERIOD, sequence_length=SEQUENCE_LENGTH_BY_PERIOD, hotspot_level_by_grid=GRID_HOTSPOT_LEVEL_BY_PERIOD, hotspot_matrix=HOTSPOT_MATRIX_BY_PERIOD, entropy_data=ENTROPY_DATA_BY_PERIOD, variance_data=VARIANCE_DATA_BY_PERIOD
    output: jsd=SIMULATION_JSD_RESULT_FILE, entropy=SIMULATION_ENTROPY_RESULT_FILE, variance=SIMULATION_VARIANCE_RESULT_FILE
    params: hotspot_level=HOTSPOT_LEVEL, max_entropy=MAX_ENTROPY, max_variance=MAX_VARIANCE, n_entropy_bin=N_ENTROPY_BIN, n_variance_bin=N_VARIANCE_BIN, repetition=REPTITION
    script: "workflow/scripts/simulation.py"

    
    
    