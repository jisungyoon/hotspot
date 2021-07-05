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


def calculate_locational_variance(given_sequences, d):
    var_array = []
    for row in given_sequences:
        counter_dict = dict(Counter(row))
        max_val = max(counter_dict.values())
        max_keys = [k for k, v in counter_dict.items() if v == max_val]
        home = max_keys[np.random.choice(len(max_keys))]
        var = 0
        for k, v in counter_dict.items():
            if k != home:
                var += (d[home][k] ** 2) * v
        var /= len(row)
        var_array.append(np.sqrt(var))

    return var_array


def get_pdf(points, bins):
    temp_counter = Counter(np.digitize(points, bins))
    temp_sum = sum(temp_counter.values())
    temp_ys = [temp_counter.get(i, 0) / temp_sum for i in range(1, len(bins))]

    return temp_ys


def generate_sequence(
    p,
    gamma,
    k,
    generated_homes,
    sequence_length,
    hotspot_level_to_grid,
    d,
):
    distribution = np.array([1 / np.power(i + 1, k) for i in range(10)])
    distribution /= sum(distribution)

    generated_sequences = []
    for home, length in zip(generated_homes, sequence_length):
        next_grids = [home]
        for i in range(length - 1):
            if np.random.random() < 1 - np.power(p, 1 + i):
                next_grid = np.random.choice(next_grids)
            else:
                level = np.random.choice(np.arange(10), p=distribution) + 1
                c = Counter(next_grids)
                m = max(c.values())
                home = np.random.choice([key for key in c if c[key] == m])
                next_grid = get_next_grid(home, gamma, d, hotspot_level_to_grid[level])
            next_grids.append(next_grid)
        generated_sequences.append(next_grids)

    return generated_sequences


def get_next_grid(home, gamma, d, target_grids):
    target_grids = [spot for spot in target_grids if spot != home]
    weight = np.array([1 / np.power(d[home][spot], gamma) for spot in target_grids])
    weight = weight / sum(weight)

    return target_grids[np.random.choice(np.arange(len(target_grids)), p=weight)]
