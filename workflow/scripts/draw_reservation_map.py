import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib
from matplotlib import cm, figure, font_manager
from pyproj import Proj, transform
from shapely.geometry import Polygon
from matplotlib.collections import PatchCollection
import s2geometry as s2

INPUT_HOTSPOT_LEVEL_BY_GRID = snakemake.input.hotspot_level_by_grid

INPUT_SHAPE_FILE = snakemake.input.shape_file
INPUT_FONT_FILE = snakemake.input.font_file

OUTPUT_HOTSPOT_MAP = snakemake.output.hotspot_map


# Import data
hotspot_info = pd.read_pickle(INPUT_HOTSPOT_LEVEL_BY_GRID)
cell_ids = list(map(int, list(hotspot_info.keys())))

# Load shape file
seoul_shp = gpd.read_file(INPUT_SHAPE_FILE)


def transform_polygon(x, mode="geo"):
    proj_UTMK = Proj(init="epsg:5178")
    proj_WGS84 = Proj(init="epsg:4326")

    utm_k_x, utm_k_y = np.array(x["geometry"].exterior.coords.xy)
    wgs_x, wgs_y = np.array(transform(proj_UTMK, proj_WGS84, utm_k_x, utm_k_y))
    if mode == "geo":
        return Polygon(zip(wgs_x, wgs_y))
    elif mode == "grid":
        new_x = (wgs_x - min_lon) / (max_lon - min_lon) * n_lon_bin
        new_y = (wgs_y - min_lat) / (max_lat - min_lat) * n_lat_bin
        return Polygon(zip(new_x, new_y))


# transform geo-encoding
seoul_shp["geometry"] = seoul_shp.apply(
    lambda x: transform_polygon(x, mode="geo"), axis=1
)

from matplotlib.patches import Polygon
patches = []
for id_ in cell_ids:
    cell = s2.S2Cell(s2.S2CellId(id_))
    vertices = []
    for i in range(0, 4):
        vertex = cell.GetVertex(i)
        latlng = s2.S2LatLng(vertex)
        vertices.append([latlng.lng().degrees(),
                         latlng.lat().degrees()])
        print(vertices)
    patches.append(Polygon(vertices, True))
    
    
prop = font_manager.FontProperties(fname=INPUT_FONT_FILE, size=22)
small_prop = font_manager.FontProperties(fname=INPUT_FONT_FILE, size=18)

f = plt.figure(figsize=(6.2, 5.6))
ax = f.add_axes([0.17, 0, 0.7, 0.7])
cax = f.add_axes([0.88, 0.04, 0.03, 0.62])

cmap = cm.get_cmap("Blues_r", 10)
seoul_shp.plot(ax=ax, color="black", edgecolor='darkgrey', alpha=0.1
)
p = PatchCollection(patches, cmap=cmap)

colors = [hotspot_info[id_] for id_ in cell_ids]
p.set_array(np.array(colors))
ax.add_collection(p)

ax.set_xlim(126.82, 127.15)
ax.set_ylim(37.44,37.70)
ax.axis("off")

norm = matplotlib.colors.BoundaryNorm(np.array(list(range(0,11))) + 0.5, 11)
cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap,
                                norm=norm,
                                ticks=[np.arange(1, 11)],
                                spacing='proportional',
                                orientation='vertical')
cbar.set_label('Hotspot Level', fontproperties=prop)
cbar.ax.invert_yaxis()
for label in cbar.ax.get_yticklabels():
    label.set_fontproperties(small_prop)

plt.savefig(OUTPUT_HOTSPOT_MAP, bbox_inches="tight")