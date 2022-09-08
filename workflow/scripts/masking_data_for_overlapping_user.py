import numpy as np

INPUT_SEQUENCE_LENGTH = snakemake.input.sequence_length
INPUT_ENTROPY = snakemake.input.entropy_data
INPUT_VARIANCE = snakemake.input.variance_data
MASKING_FILE = snakemake.input.masking

OUTPUT_SEQUENCE_LENGTH = snakemake.output.masked_sequence_length
OUTPUT_ENTROPY = snakemake.output.masked_entropy_data
OUTPUT_VARIANCE = snakemake.output.masked_variance_data


mask = np.load(MASKING_FILE)

sequence_length = np.load(INPUT_SEQUENCE_LENGTH)
entropy = np.load(INPUT_ENTROPY)
variance = np.load(INPUT_VARIANCE)

np.save(OUTPUT_SEQUENCE_LENGTH, sequence_length[mask])
np.save(OUTPUT_ENTROPY, entropy[mask])
np.save(OUTPUT_VARIANCE, variance[mask])
