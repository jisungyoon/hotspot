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


def generate_sequence(p, k, generated_homes, sequence_length, grid_to_hotspot_level, hotspot_level_to_grid, transition_matrix):
    generated_sequences = []
    for home, length in zip(generated_homes, sequence_length):
        home_hotspot_level = int(grid_to_hotspot_level[home])
        temp_trans_prob = transition_matrix[home_hotspot_level - 1]
        next_hot_levels = np.random.choice(np.arange(10), size=int(length/2), replace=True, p=temp_trans_prob) + 1
        next_grids = []
        
        for level in next_hot_levels:
            if home_hotspot_level == level:
                if np.random.random() < p:
                    next_grid = home
                else:
                    next_grid = get_next_grid(home, k, hotspot_level_to_grid[level])
            else:
                next_grid = get_next_grid(home, k, hotspot_level_to_grid[level])
            next_grids.append(next_grid)
        
        temp_sequnece = []
        for i in range(length):
            if i % 2 == 0:
                temp_sequnece.append(home)
            else:
                temp_sequnece.append(next_grids.pop())
        generated_sequences.append(temp_sequnece)
        
    return generated_sequences
        
def get_next_grid(home, k, target_grids):
    target_grids = [spot for spot in target_grids if spot != home]
    weight = np.array([1/np.power(np.sqrt((spot[0] - home[0]) ** 2 + (spot[1] - home[1]) ** 2), k) for spot in target_grids])
    weight = weight/sum(weight)
    
    return target_grids[np.random.choice(np.arange(len(target_grids)), p=weight)]