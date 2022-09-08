from os.path import join as j
import numpy as np

configfile: "workflow/config.yaml"

# Define directory and configuration for the project
DATA_DIR = config["data_dir"]
RAW_DIR = j(DATA_DIR, "raw")
DERIVED_DIR = j(DATA_DIR, "derived")
SHAPE_FILE_DIR = j(DATA_DIR, "shape_file", "seoul")
SIMULATION_DIR = config["simulation_dir"]
FIGURE_DIR = config["figs_dir"]
ASSETS_DIR = config["assets_dir"] 
FIG_FORMAT = 'pdf'

###############################################################################
# RAW_DATA
###############################################################################
PERIODS = ["before", "after"]
RAW_SEQUENCE_DATA = j(RAW_DIR, "{period}_data.csv")
META_DATA = j(RAW_DIR, "meta_data.csv")
SHAPE_FILE = j(SHAPE_FILE_DIR, "seoul.shp")

###############################################################################
# DERIVED_DATA
###############################################################################
DERIVED_DATA_BY_PERIOD = j(DERIVED_DIR, "{period}")
MASKING_FOR_OVERLAPPING_USER = j(DERIVED_DATA_BY_PERIOD, "mask.npy")
SEQUENCE_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "sequence.npy")

GRID_RESERVATION_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "reservation_by_grid.pkl")
GRID_HOTSPOT_LEVEL_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "hotspot_level_by_grid.pkl")

VENDOR_TO_GRID_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "vendor_to_grid.pkl")
VENDOR_TO_HOTSPOT_LEVEL_BY_PERIOD = j(
    DERIVED_DATA_BY_PERIOD, "vendor_to_hotspot_level.pkl"
)

PSI_RESULT_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "psi.csv")
OVERLAP_PSI_RESULT_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "overlap_psi.csv")

DISTANCE = j(DERIVED_DATA_BY_PERIOD, "distance.pkl")
HOME_PDF_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "home_pdf.pkl")

SEQUENCE_LENGTH_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "sequence_length.npy")
ENTROPY_DATA_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "entropy_data.npy")
VARIANCE_DATA_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "varaince_data.npy")

# For overlapping user
MASKED_SEQUENCE_LENGTH_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "masked_sequence_length.npy")
MASKED_ENTROPY_DATA_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "masked_entropy_data.npy")
MASKED_VARIANCE_DATA_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "masked_varaince_data.npy")

PS_DATA_BY_PERIOD =j(DERIVED_DATA_BY_PERIOD, "p.npy")

###############################################################################
# SIMULATION_RESULT_FOR_THE_MAIN_RESULT
###############################################################################
OVERLAP_SIMULATION_RESULT_DIR = j(SIMULATION_DIR, "{period}")

SIMULATION_RESULT_DIR_BY_PK = j(OVERLAP_SIMULATION_RESULT_DIR, "find_pk_multi")

SIMULATION_JSD_RESULT_FILE_BY_PK = j(
    SIMULATION_RESULT_DIR_BY_PK, "result", "result_{p}_{k}.csv"
)
SIMULATION_ENTROPY_RESULT_FILE_BY_PK = j(
    SIMULATION_RESULT_DIR_BY_PK, "entropy", "entropy_{p}_{k}.npy"
)
SIMULATION_VARIANCE_RESULT_FILE_BY_PK = j(
    SIMULATION_RESULT_DIR_BY_PK, "variance", "variance_{p}_{k}.npy"
)

SIMULATION_RESULT_DIR_BY_GAMMA = j(OVERLAP_SIMULATION_RESULT_DIR, "find_gamma_final")

SIMULATION_JSD_RESULT_FILE_BY_GAMMA = j(
    SIMULATION_RESULT_DIR_BY_GAMMA, "result", "result_{gamma}.csv"
)
SIMULATION_ENTROPY_RESULT_FILE_BY_GAMMA = j(
    SIMULATION_RESULT_DIR_BY_GAMMA, "entropy", "entropy_{gamma}.npy"
)
SIMULATION_VARIANCE_RESULT_FILE_BY_GAMMA = j(
    SIMULATION_RESULT_DIR_BY_GAMMA, "variance", "variance_{gamma}.npy"
)

###############################################################################
# ASSETS
###############################################################################
NORMAL_FONT_PATH = j(ASSETS_DIR, "Helvetica.ttf")

###############################################################################
# FIGS
###############################################################################
FIGS_BY_PERIOD = j(FIGURE_DIR, "{period}")
HOTSPOT_MAP_BY_PERIOD = j(FIGS_BY_PERIOD, "hotspot_map.{}".format(FIG_FORMAT))

HOTSPOT_MATRIX_FIG_BY_PERIOD = j(FIGS_BY_PERIOD, "hotspot_matrix.{}".format(FIG_FORMAT))
NULL_HOTSPOT_MATRIX_FIG_BY_PERIOD = j(FIGS_BY_PERIOD, "null_hotspot_matrix.{}".format(FIG_FORMAT))
RATIO_FIG_BY_PERIOD = j(FIGS_BY_PERIOD, "ratio_matrix.{}".format(FIG_FORMAT))

HOME_VALIDATION_BY_PERIOD = j(FIGS_BY_PERIOD, "home_validation.{}".format(FIG_FORMAT))

ENT_VAR_BY_PERIOD = j(FIGS_BY_PERIOD, "entropy_variance_data.{}".format(FIG_FORMAT))
OVERLAP_ENT_VAR = j(FIGURE_DIR, "entropy_variance_data.{}".format(FIG_FORMAT))

SIMULATION_FIG_BY_PK = j(FIGURE_DIR, "simulation_by_pk_overlap_multi.{}".format(FIG_FORMAT))
SIMULATION_FIG_BY_GAMMA = j(FIGURE_DIR, "simulation_by_gamma_overlap_multi.{}".format(FIG_FORMAT))

###############################################################################
# PARAMS FOR SEOUL
###############################################################################
N_LAT_BIN = 30
N_LON_BIN = 30

MAX_LAT = 37.69
MIN_LAT = 37.45
MAX_LON = 127.16
MIN_LON = 126.80

HOTSPOT_LEVEL = 10

MAX_ENTROPY = 2.1
MAX_VARIANCE = 18

N_ENTROPY_BIN = 30
N_VARIANCE_BIN = 30

###############################################################################
# PARAMS FOR SIMULATION
###############################################################################
ps = [np.round(x, 3) for x in np.linspace(0, 1, 51)]
ks = [np.round(x, 3) for x in np.linspace(0, 3, 121)]
gammas = [np.round(x, 3) for x in np.linspace(0, 5, 201)]
REPETITION = 10

###############################################################################
# ANALYSIS PIPELINE
###############################################################################
rule all:
    input:
        expand(HOTSPOT_MAP_BY_PERIOD, period=PERIODS),
        expand(HOTSPOT_MATRIX_FIG_BY_PERIOD, period=PERIODS),
        expand(NULL_HOTSPOT_MATRIX_FIG_BY_PERIOD, period=PERIODS),
        expand(RATIO_FIG_BY_PERIOD, period=PERIODS),
        OVERLAP_ENT_VAR,
        SIMULATION_FIG_BY_PK,
        SIMULATION_FIG_BY_GAMMA
              

rule generate_sequence:
    input:
        RAW_SEQUENCE_DATA,
    output:
        sequence=SEQUENCE_BY_PERIOD,
        sequence_length=SEQUENCE_LENGTH_BY_PERIOD,
    script:
        "workflow/scripts/generate_sequence.py"


rule hotspot_analysis:
    input:
        raw=RAW_SEQUENCE_DATA,
        meta=META_DATA,
    output:
        reservation_by_grid=GRID_RESERVATION_BY_PERIOD,
        hotspot_level_by_grid=GRID_HOTSPOT_LEVEL_BY_PERIOD,
        vendor_to_grid=VENDOR_TO_GRID_BY_PERIOD,
        vendor_to_hotspot_level=VENDOR_TO_HOTSPOT_LEVEL_BY_PERIOD,
    params:
        n_lat_bin=N_LAT_BIN,
        n_lon_bin=N_LON_BIN,
        max_lat=MAX_LAT,
        min_lat=MIN_LAT,
        max_lon=MAX_LON,
        min_lon=MIN_LON,
        hotspot_level=HOTSPOT_LEVEL,
    script:
        "workflow/scripts/make_grids_and_hotspot_anlaysis.py"
        
  
rule draw_hotspot_map:
    input:
        hotspot_level_by_grid=GRID_HOTSPOT_LEVEL_BY_PERIOD,
        shape_file=SHAPE_FILE,
        font_file=NORMAL_FONT_PATH,
    output:
        hotspot_map=HOTSPOT_MAP_BY_PERIOD
    script:
        "workflow/scripts/draw_reservation_map.py"      
        

rule calculate_transition_matrix_and_psi_overlap:
    input:
        sequence=SEQUENCE_BY_PERIOD,
        masking=MASKING_FOR_OVERLAPPING_USER,
        hotspot_level_by_grid=GRID_HOTSPOT_LEVEL_BY_PERIOD,
        vendor_to_grid=VENDOR_TO_GRID_BY_PERIOD,
        font_file=NORMAL_FONT_PATH,
    output:
        hotspot_matrix_fig=HOTSPOT_MATRIX_FIG_BY_PERIOD,
        null_hotspot_matrix_fig=NULL_HOTSPOT_MATRIX_FIG_BY_PERIOD,
        ratio_fig=RATIO_FIG_BY_PERIOD,
        psi_result=PSI_RESULT_BY_PERIOD,
    params:
        hotspot_level=HOTSPOT_LEVEL,
    script:
        "workflow/scripts/calculate_transition_matrix_and_psi_overlap.py"


rule calculate_distance_and_travel_time:
    input:
        hotspot_data=GRID_HOTSPOT_LEVEL_BY_PERIOD,
    output:
        distance=DISTANCE,
    script:
        "workflow/scripts/calculate_geodistance.py"


rule calculate_and_validate_home_distribution:
    input:
        sequence=SEQUENCE_BY_PERIOD,
        vendor_to_grid=VENDOR_TO_GRID_BY_PERIOD,
        reservation_by_grid=GRID_RESERVATION_BY_PERIOD,
        font_file=NORMAL_FONT_PATH,
    output:
        home_pdf=HOME_PDF_BY_PERIOD,
        home_validation_fig=HOME_VALIDATION_BY_PERIOD,
    script:
        "workflow/scripts/calculate_and_validate_home_distribution.py"


rule calculate_entropy_variance_data:
    input:
        sequence=SEQUENCE_BY_PERIOD,
        vendor_to_grid=VENDOR_TO_GRID_BY_PERIOD,
        grid_to_hotspot_level=GRID_HOTSPOT_LEVEL_BY_PERIOD,
        font_file=NORMAL_FONT_PATH,
        dist=DISTANCE,
    output:
        entropy=ENTROPY_DATA_BY_PERIOD,
        variance=VARIANCE_DATA_BY_PERIOD,
        ent_var_fig=ENT_VAR_BY_PERIOD,
    params:
        hotspot_level=HOTSPOT_LEVEL,
        max_entropy=MAX_ENTROPY,
        max_variance=MAX_VARIANCE,
        n_entropy_bin=N_ENTROPY_BIN,
        n_variance_bin=N_VARIANCE_BIN,
    script:
        "workflow/scripts/calculate_entropy_and_variance_data.py"


rule masking_data_for_overlapping_user:
    input:
        sequence_length=SEQUENCE_LENGTH_BY_PERIOD,
        entropy_data=ENTROPY_DATA_BY_PERIOD,
        variance_data=VARIANCE_DATA_BY_PERIOD,
        masking=MASKING_FOR_OVERLAPPING_USER,
    output:
        masked_sequence_length=MASKED_SEQUENCE_LENGTH_BY_PERIOD,
        masked_entropy_data=MASKED_ENTROPY_DATA_BY_PERIOD,
        masked_variance_data=MASKED_VARIANCE_DATA_BY_PERIOD
    script:
        "workflow/scripts/masking_data_for_overlapping_user.py"


rule draw_overlap_ent_var_fig:
    input:
        entropy=expand(MASKED_ENTROPY_DATA_BY_PERIOD, period=PERIODS),
        variance=expand(MASKED_VARIANCE_DATA_BY_PERIOD, period=PERIODS),
        font_file=NORMAL_FONT_PATH,
    output:
        overlap_ent_var_fig=OVERLAP_ENT_VAR,
    params:
        hotspot_level=HOTSPOT_LEVEL,
        max_entropy=MAX_ENTROPY,
        max_variance=MAX_VARIANCE,
        n_entropy_bin=N_ENTROPY_BIN,
        n_variance_bin=N_VARIANCE_BIN,
    script:
        "workflow/scripts/draw_overlapped_ent_var.py"
             
    
###############################################################################
# SIMULATION
###############################################################################
rule calculate_p:
    input:
        sequence=SEQUENCE_BY_PERIOD,
        masking=MASKING_FOR_OVERLAPPING_USER,
        vendor_to_grid=VENDOR_TO_GRID_BY_PERIOD,
        grid_to_hotspot_level=GRID_HOTSPOT_LEVEL_BY_PERIOD,
    output:
        p=PS_DATA_BY_PERIOD
    params:
        hotspot_level=HOTSPOT_LEVEL,
    script:
        "workflow/scripts/calculate_p.py"
        
        
rule simulation_by_p_and_k:
    input:
        home_pdf=HOME_PDF_BY_PERIOD,
        sequence_length=MASKED_SEQUENCE_LENGTH_BY_PERIOD,
        hotspot_level_by_grid=GRID_HOTSPOT_LEVEL_BY_PERIOD,
        entropy_data=MASKED_ENTROPY_DATA_BY_PERIOD,
        variance_data=MASKED_VARIANCE_DATA_BY_PERIOD,
        p_data=PS_DATA_BY_PERIOD,
        dist=DISTANCE,
    output:
        jsd=SIMULATION_JSD_RESULT_FILE_BY_PK,
        entropy=SIMULATION_ENTROPY_RESULT_FILE_BY_PK,
        variance=SIMULATION_VARIANCE_RESULT_FILE_BY_PK,
    params:
        hotspot_level=HOTSPOT_LEVEL,
        max_entropy=MAX_ENTROPY,
        max_variance=MAX_VARIANCE,
        n_entropy_bin=N_ENTROPY_BIN,
        n_variance_bin=N_VARIANCE_BIN,
        repetition=10,
    script:
        "workflow/scripts/simulation.py"
 
 
rule draw_simulation_result_by_p_and_k:
    input:
        before_files=expand(
            SIMULATION_JSD_RESULT_FILE_BY_PK, period=["before"], p=ps, k=ks
        ),
        after_files=expand(
            SIMULATION_JSD_RESULT_FILE_BY_PK, period=["after"], p=ps, k=ks
        ),
        font_file=NORMAL_FONT_PATH,
    output:
        simulation_fig_by_pk=SIMULATION_FIG_BY_PK,
    params:
        ps=ps,
        ks=ks,
    script:
        "workflow/scripts/draw_simulation_result_hotspot_entropy.py"
        
        
rule simulation_by_gamma:
    input:
        home_pdf=HOME_PDF_BY_PERIOD,
        sequence_length=MASKED_SEQUENCE_LENGTH_BY_PERIOD,
        hotspot_level_by_grid=GRID_HOTSPOT_LEVEL_BY_PERIOD,
        entropy_data=MASKED_ENTROPY_DATA_BY_PERIOD,
        variance_data=MASKED_VARIANCE_DATA_BY_PERIOD,
        p_data=PS_DATA_BY_PERIOD,
        dist=DISTANCE,
    output:
        jsd=SIMULATION_JSD_RESULT_FILE_BY_GAMMA,
        entropy=SIMULATION_ENTROPY_RESULT_FILE_BY_GAMMA,
        variance=SIMULATION_VARIANCE_RESULT_FILE_BY_GAMMA,
    params:
        hotspot_level=HOTSPOT_LEVEL,
        max_entropy=MAX_ENTROPY,
        max_variance=MAX_VARIANCE,
        n_entropy_bin=N_ENTROPY_BIN,
        n_variance_bin=N_VARIANCE_BIN,
        repetition=REPETITION,
    script:
        "workflow/scripts/simulation_find_gamma.py"
              
        
rule draw_simulation_result_by_gamma:
    input:
        before_files=expand(
            SIMULATION_JSD_RESULT_FILE_BY_GAMMA, period=["before"], gamma=gammas
        ),
        after_files=expand(
            SIMULATION_JSD_RESULT_FILE_BY_GAMMA, period=["after"], gamma=gammas
        ),
        font_file=NORMAL_FONT_PATH,
    output:
        simulation_fig_by_gamma=SIMULATION_FIG_BY_GAMMA,
    params:
        gammas=gammas,
    script:
        "workflow/scripts/draw_simulation_result_locational_variance.py"