import datetime
import time
from collections import defaultdict

import numpy as np
import pandas as pd
import s2geometry as s2
from geopy.distance import great_circle
from tqdm import tqdm

HOTSPOT_DATA = snakemake.input.hotspot_data

OUTPUT_DISTANCE = snakemake.output.distance


cell_ids = list(map(int, list(pd.read_pickle(HOTSPOT_DATA).keys())))


def compute_geo_distance(u, v):
    """Compute geographic distance (great circle), between two
       sets of coordinates. Outputs NaN when any distnace is NaN
    arguments:
    u -- First coordinate, list in form [<lat>, <lon>]
    v -- Second coordinate, list in form [<lat>, <lon>]
    """
    if np.isnan([u[0], u[1], v[0], v[1]]).any():
        return np.nan
    else:
        return great_circle(u, v).kilometers


def default_to_regular(d):
    if isinstance(d, defaultdict):
        d = {k: default_to_regular(v) for k, v in d.items()}
    return d


center = {}
for id_ in cell_ids:
    cell = s2.S2Cell(s2.S2CellId(id_))
    vertices = []
    for i in range(0, 4):
        vertex = cell.GetVertex(i)
        latlng = s2.S2LatLng(vertex)
        vertices.append([latlng.lng().degrees(), latlng.lat().degrees()])
    center[id_] = np.mean(vertices, axis=0)

cell_ks = list(center.keys())
distance = defaultdict(lambda: defaultdict(int))

for origin in tqdm(cell_ks):
    origin_cord = (center[origin][1], center[origin][0])
    for dest in cell_ks:
        if origin != dest:
            dest_cord = (center[dest][1], center[dest][0])
            distance[origin][dest] = compute_geo_distance(origin_cord, dest_cord)

pd.to_pickle(default_to_regular(distance), OUTPUT_DISTANCE)
