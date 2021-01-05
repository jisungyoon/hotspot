import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm, figure, font_manager
from matplotlib.colors import LogNorm
from pyproj import Proj, transform
from shapely.geometry import Polygon

INPUT_RESERVATION_BY_GRID = snakemake.input.reservation_by_grid
INPUT_HOTSPOT_LEVEL_BY_GRID = snakemake.input.hotspot_level_by_grid

INPUT_SHAPE_FILE = snakemake.input.shape_file
INPUT_FONT_FILE = snakemake.input.font_file

OUTPUT_RESERVATION_MAP = snakemake.output.reservation_map
OUTPUT_HOTSPOT_MAP = snakemake.output.hotspot_map


n_lat_bin = snakemake.params.n_lat_bin
n_lon_bin = snakemake.params.n_lon_bin

max_lat = snakemake.params.max_lat
min_lat = snakemake.params.min_lat
max_lon = snakemake.params.max_lon
min_lon = snakemake.params.min_lon


# Import data
reservation_by_grid = np.load(INPUT_RESERVATION_BY_GRID)
hotspot_level_by_grid = np.load(INPUT_HOTSPOT_LEVEL_BY_GRID)

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
    lambda x: transform_polygon(x, mode="grid"), axis=1
)

# draw reservation map
prop = font_manager.FontProperties(fname=INPUT_FONT_FILE, size=22)
small_prop = font_manager.FontProperties(fname=INPUT_FONT_FILE, size=18)

f = plt.figure(figsize=(6.2, 5.6))
ax = f.add_axes([0.17, 0.02, 0.72, 0.79])
axcolor = f.add_axes([0.93, 0.02, 0.03, 0.79])

seoul_shp.plot(ax=ax, color="black", edgecolor="grey", alpha=0.1)
im = ax.pcolormesh(reservation_by_grid.T, cmap=cm.Blues, norm=LogNorm())
ax.set_xlim([0, 30])
ax.set_ylim([0, 30])
ax.axis("off")

cbar = f.colorbar(im, cax=axcolor)
for label in cbar.ax.get_yticklabels():
    label.set_fontproperties(small_prop)
cbar.ax.set_ylabel("Reservation", fontproperties=prop)

plt.savefig(OUTPUT_RESERVATION_MAP, bbox_inches="tight")


# draw hotspot map
masked_hotspot_level_by_grid = np.ma.masked_where(
    hotspot_level_by_grid == 11, hotspot_level_by_grid
)

f = plt.figure(figsize=(6.2, 5.6))
ax = f.add_axes([0.17, 0.02, 0.72, 0.79])
axcolor = f.add_axes([0.90, 0.02, 0.03, 0.79])

cmap = cm.get_cmap("Blues_r", 10)
seoul_shp.plot(ax=ax, color="black", edgecolor="grey", alpha=0.1)
im = ax.pcolormesh(masked_hotspot_level_by_grid.T, cmap=cmap, vmin=0.5, vmax=10.5)
ax.set_xlim([0, 30])
ax.set_ylim([0, 30])
ax.axis("off")

cbar = f.colorbar(im, cax=axcolor, ticks=np.arange(1, 11))
cbar.ax.invert_yaxis()
for label in cbar.ax.get_yticklabels():
    label.set_fontproperties(small_prop)
cbar.ax.set_ylabel("Hotspot Level", fontproperties=prop)

plt.savefig(OUTPUT_HOTSPOT_MAP, bbox_inches="tight")

f.show()
