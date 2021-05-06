from collections import Counter

import numpy as np
import pandas as pd

INPUT = snakemake.input[0]
OUTPUT_SEQUENCE = snakemake.output.sequence
OUTPUT_SEQUENCE_LENGTH = snakemake.output.sequence_length


data = pd.read_csv(INPUT, sep="|", parse_dates=["dt_reservation", "dt_entrance"])
data = data[data.order_cnt > 0]  # using only valid reservation
data = data[data.uno2 != 0.0]  # remove Non-Members

sequences = []
for user_id, row in data.groupby("uno2"):
    if len(row) > 1:
        sorted_row = row.sort_values("dt_entrance")
        sequences.append(list(sorted_row["yg_vendor_code"]))
sequence_length = list(map(len, sequences))

np.save(OUTPUT_SEQUENCE, sequences)
np.save(OUTPUT_SEQUENCE_LENGTH, sequence_length)
