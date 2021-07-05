from os.path import join as j
import numpy as np


configfile: "workflow/config.yaml"


API_KEY = config["api_key"]
DATA_DIR = config["data_dir"]
SIMULATION_DIR = config["simulation_dir"]
RAW_DIR = j(DATA_DIR, "raw")
DERIVED_DIR = j(DATA_DIR, "derived")
SHAPE_FILE_DIR = j(DATA_DIR, "shape_file", "seoul")

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
SEQUENCE_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "sequence.npy")
SEQUENCE_LENGTH_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "sequence_length.npy")

GRID_RESERVATION_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "reservation_by_grid.pkl")
GRID_HOTSPOT_LEVEL_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "hotspot_level_by_grid.pkl")

VENDOR_TO_GRID_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "vendor_to_grid.pkl")
VENDOR_TO_HOTSPOT_LEVEL_BY_PERIOD = j(
    DERIVED_DATA_BY_PERIOD, "vendor_to_hotspot_level.pkl"
)

PSI_RESULT_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "psi.csv")

DISTANCE = j(DERIVED_DATA_BY_PERIOD, "distance.pkl")
HOME_PDF_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "home_pdf.pkl")

ENTROPY_DATA_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "entropy_data.npy")
VARIANCE_DATA_BY_PERIOD = j(DERIVED_DATA_BY_PERIOD, "varaince_data.npy")

###############################################################################
# SIMULATION_RESULT
###############################################################################
SIMULATION_RESULT_DIR = j(SIMULATION_DIR, "{period}")

SIMULATION_RESULT_DIR_BY_PK = j(SIMULATION_RESULT_DIR, "find_pk")

SIMULATION_JSD_RESULT_FILE_BY_PK = j(
    SIMULATION_RESULT_DIR_BY_PK, "result", "result_{p}_{k}.csv"
)
SIMULATION_ENTROPY_RESULT_FILE_BY_PK = j(
    SIMULATION_RESULT_DIR_BY_PK, "entropy", "entropy_{p}_{k}.npy"
)
SIMULATION_VARIANCE_RESULT_FILE_BY_PK = j(
    SIMULATION_RESULT_DIR_BY_PK, "variance", "variance_{p}_{k}.npy"
)


SIMULATION_RESULT_DIR_BY_GAMMA = j(SIMULATION_RESULT_DIR, "find_gamma")

SIMULATION_JSD_RESULT_FILE_BY_GAMMA = j(
    SIMULATION_RESULT_DIR_BY_GAMMA, "result", "result_{gamma}.csv"
)
SIMULATION_ENTROPY_RESULT_FILE_BY_GAMMA = j(
    SIMULATION_RESULT_DIR_BY_GAMMA, "entropy", "entropy_{gamma}.npy"
)
SIMULATION_VARIANCE_RESULT_FILE_BY_GAMMA = j(
    SIMULATION_RESULT_DIR_BY_GAMMA, "variance", "variance_{gamma}.npy"
)

SIMULATION_RESULT_DIR_BY_GAMMA_WITHOUT_POPULARITY = j(SIMULATION_RESULT_DIR, "without_popularity")

SIMULATION_JSD_RESULT_FILE_BY_GAMMA_WITHOUT_POPULARITY = j(
    SIMULATION_RESULT_DIR_BY_GAMMA_WITHOUT_POPULARITY, "result", "result_{gamma}.csv"
)
SIMULATION_ENTROPY_RESULT_FILE_BY_GAMMA_WITHOUT_POPULARITY = j(
    SIMULATION_RESULT_DIR_BY_GAMMA_WITHOUT_POPULARITY, "entropy", "entropy_{gamma}.npy"
)
SIMULATION_VARIANCE_RESULT_FILE_BY_GAMMA_WITHOUT_POPULARITY = j(
    SIMULATION_RESULT_DIR_BY_GAMMA_WITHOUT_POPULARITY, "variance", "variance_{gamma}.npy"
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
OVERLAP_ENT_VAR = j(FIGURE_DIR, "entropy_variance_data_overlap.{}".format(FIG_FORMAT))

SIMULATION_FIG_BY_PK = j(FIGURE_DIR, "simulation_by_pk.{}".format(FIG_FORMAT))
SIMULATION_FIG_BY_GAMMA = j(FIGURE_DIR, "simulation_by_gamma.{}".format(FIG_FORMAT))
SIMULATION_FIG_BY_GAMMA_WITHOUT_POPULARITY = j(FIGURE_DIR, "simulation_by_gamma_without_popularity.{}".format(FIG_FORMAT))


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

MAX_ENTROPY = 2.1
MAX_VARIANCE = 18

N_ENTROPY_BIN = 30
N_VARIANCE_BIN = 30

###############################################################################
# PARAMS for simulation
###############################################################################
ps = [np.round(x, 3) for x in np.linspace(0, 1, 51)]
ks = [np.round(x, 3) for x in np.linspace(0, 3, 121)]

REPTITION = 10
gammas = [np.round(x, 2) for x in np.linspace(0, 5, 201)]


rule all:
    input:
        expand(HOTSPOT_MAP_BY_PERIOD, period=PERIODS),
        expand(HOTSPOT_MATRIX_FIG_BY_PERIOD, period=PERIODS),
        expand(NULL_HOTSPOT_MATRIX_FIG_BY_PERIOD, period=PERIODS),
        OVERLAP_ENT_VAR,
        SIMULATION_FIG_BY_PK,
        SIMULATION_FIG_BY_GAMMA,
        SIMULATION_FIG_BY_GAMMA_WITHOUT_POPULARITY


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
        


rule calculate_transition_matrix_and_psi:
    input:
        sequence=SEQUENCE_BY_PERIOD,
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
        "workflow/scripts/calculate_transition_matrix_and_psi.py"


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


rule draw_overlap_ent_var_fig:
    input:
        entropy=expand(ENTROPY_DATA_BY_PERIOD, period=PERIODS),
        variance=expand(VARIANCE_DATA_BY_PERIOD, period=PERIODS),
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


rule simulation_by_p_and_k:
    input:
        home_pdf=HOME_PDF_BY_PERIOD,
        sequence_length=SEQUENCE_LENGTH_BY_PERIOD,
        hotspot_level_by_grid=GRID_HOTSPOT_LEVEL_BY_PERIOD,
        entropy_data=ENTROPY_DATA_BY_PERIOD,
        variance_data=VARIANCE_DATA_BY_PERIOD,
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
        repetition=REPTITION,
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
        sequence_length=SEQUENCE_LENGTH_BY_PERIOD,
        hotspot_level_by_grid=GRID_HOTSPOT_LEVEL_BY_PERIOD,
        entropy_data=ENTROPY_DATA_BY_PERIOD,
        variance_data=VARIANCE_DATA_BY_PERIOD,
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
        repetition=REPTITION,
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
        
        
rule simulation_by_gamma_without_popularity:
    input:
        home_pdf=HOME_PDF_BY_PERIOD,
        sequence_length=SEQUENCE_LENGTH_BY_PERIOD,
        hotspot_level_by_grid=GRID_HOTSPOT_LEVEL_BY_PERIOD,
        entropy_data=ENTROPY_DATA_BY_PERIOD,
        variance_data=VARIANCE_DATA_BY_PERIOD,
        dist=DISTANCE,
    output:
        jsd=SIMULATION_JSD_RESULT_FILE_BY_GAMMA_WITHOUT_POPULARITY,
        entropy=SIMULATION_ENTROPY_RESULT_FILE_BY_GAMMA_WITHOUT_POPULARITY,
        variance=SIMULATION_VARIANCE_RESULT_FILE_BY_GAMMA_WITHOUT_POPULARITY,
    params:
        hotspot_level=HOTSPOT_LEVEL,
        max_entropy=MAX_ENTROPY,
        max_variance=MAX_VARIANCE,
        n_entropy_bin=N_ENTROPY_BIN,
        n_variance_bin=N_VARIANCE_BIN,
        repetition=REPTITION,
    script:
        "workflow/scripts/simulation_find_gamma_without_popularity.py"
        
rule draw_simulation_result_by_gamma_without_popularity:
    input:
        before_files=expand(
            SIMULATION_JSD_RESULT_FILE_BY_GAMMA_WITHOUT_POPULARITY, period=["before"], gamma=gammas
        ),
        after_files=expand(
            SIMULATION_JSD_RESULT_FILE_BY_GAMMA_WITHOUT_POPULARITY, period=["after"], gamma=gammas
        ),
        font_file=NORMAL_FONT_PATH,
    output:
        simulation_fig_by_gamma=SIMULATION_FIG_BY_GAMMA_WITHOUT_POPULARITY,
    params:
        gammas=gammas,
    script:
        "workflow/scripts/draw_simulation_result_locational_variance.py"
