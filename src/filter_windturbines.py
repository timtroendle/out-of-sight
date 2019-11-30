import json
from functools import partial

import rasterio
import pyproj
import shapely.geometry
import shapely.ops


WGS84 = "EPSG:4326"
EPSG_3035 = "EPSG:3035"
REPLACEMENT_VALUE = 20
RASTER_SIZE_M = 2.5
TURBINE_AREA_M = 100 # ASSUME no buildings within 50 m distance to existing turbine


def filter_turbines(path_to_esm, path_to_turbines, path_to_output):
    with rasterio.open(path_to_esm) as esm_src:
        esm = esm_src.read(1)
        meta = esm_src.meta
    with open(path_to_turbines, "r") as src:
        turbines = json.load(src)
    for turbine_loc in [(turbine["lat"], turbine["lon"]) for turbine in turbines["elements"]]:
        x, y = esm_src.index(*transform_coordinates(turbine_loc[1], turbine_loc[0], WGS84, EPSG_3035))
        shift = int(TURBINE_AREA_M / 2 / RASTER_SIZE_M)
        esm[x - shift:x + shift, y - shift:y + shift] = 0

    meta["compress"] = "lzw"
    with rasterio.open(path_to_output, 'w', **meta) as dst:
        dst.write(esm, 1)


def transform_coordinates(x, y, from_epsg, to_epsg):
    """Tranforms coordinates from one coordinate reference system to the other."""
    point = shapely.geometry.Point(x, y)

    project = partial(
        pyproj.transform,
        pyproj.Proj(init=from_epsg),
        pyproj.Proj(init=to_epsg))

    transformed_point = shapely.ops.transform(project, point)
    return transformed_point.x, transformed_point.y


if __name__ == "__main__":
    filter_turbines(
        path_to_esm=snakemake.input.esm,
        path_to_turbines=snakemake.input.turbines,
        path_to_output=snakemake.output[0]
    )
