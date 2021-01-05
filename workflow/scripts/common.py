from collections import Counter, defaultdict

import numpy as np
from scipy.stats import entropy


def calculate_entropy(given_sequences, grid_to_hotspot_level, hotspot_level=10):
    entropy_array = []
    for row in given_sequences:
        hotspot_row = list(map(grid_to_hotspot_level.get, row))
        temp_vec = np.zeros(hotspot_level)
        temp_dict = dict(Counter(hotspot_row))
        temp_sum = sum(temp_dict.values())
        for k, v in temp_dict.items():
            temp_vec[int(k - 1)] = v / temp_sum
        entropy_array.append(entropy(temp_vec))
    entropy_array = np.array(entropy_array)

    return entropy_array


def calculate_locational_variance(given_sequences):
    var_array = []
    for row in given_sequences:
        counter_dict = dict(Counter(row))
        max_val = max(counter_dict.values())
        max_keys = [k for k, v in counter_dict.items() if v == max_val]
        home = max_keys[np.random.choice(len(max_keys))]
        x_vars = 0
        y_vars = 0
        for k, v in counter_dict.items():
            x_vars += ((k[0] - home[0]) ** 2) * v
            y_vars += ((k[1] - home[1]) ** 2) * v
        var = np.sqrt((x_vars + y_vars) / len(row))
        var_array.append(var)
    var_array = np.array(var_array)

    return var_array


def get_pdf(points, bins):
    temp_counter = Counter(np.digitize(points, bins))
    temp_sum = sum(temp_counter.values())
    temp_ys = [temp_counter.get(i, 0) / temp_sum for i in range(1, len(bins))]

    return temp_ys
