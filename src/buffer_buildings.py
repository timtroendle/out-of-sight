import gc

import numpy as np
import rasterio
import rasterio.features
import shapely
import shapely.geometry


BUILDING_CODE = 50
RESOLUTION_IN_M = 2.5
SCALING_FACTOR = 10 # results in data 100 times smaller


def buffer_buildings(path_to_esm, distance, path_to_output):
    with rasterio.open(path_to_esm) as src:
        esm = src.read(1)
        meta = src.meta

    buffered = buffer_buildings_using_circles(esm, buffer_in_m=distance, transform=meta["transform"])
    del esm
    gc.collect() # immediately remove unnecessary data to avoid memory peaks

    smaller, small_meta = downsample(buffered, meta, SCALING_FACTOR)
    del buffered
    gc.collect()

    with rasterio.open(path_to_output, 'w', **small_meta) as dst:
        dst.write(smaller)


def buffer_buildings_using_circles(esm, buffer_in_m, transform):
    mask = esm == BUILDING_CODE
    shapes = rasterio.features.shapes(
        mask.astype(np.uint8),
        mask=mask,
        transform=transform
    )
    shapes = ((shapely.geometry.shape(shape).buffer(buffer_in_m), value)
              for shape, value in shapes)
    return rasterio.features.rasterize(
        shapes,
        out_shape=esm.shape,
        fill=0.0,
        transform=transform,
        dtype=np.uint8
    )


def downsample(data, meta, scale):
    meta = meta.copy()
    with rasterio.MemoryFile() as memfile:
        meta["dtype"] = np.float32
        meta["nodata"] = -1
        width = meta["width"]
        with memfile.open(**meta) as dataset:
            dataset.write(data.astype(np.float32, copy=False), 1)
            smaller = dataset.read(
                out_shape=(width // scale, width // scale),
                resampling=rasterio.enums.Resampling.average
            )
            small_meta = meta.copy()
            t = small_meta["transform"]
            small_meta["transform"] = rasterio.Affine(t.a * scale, t.b, t.c, t.d, t.e * scale, t.f)
            small_meta["height"] = small_meta["height"] // scale
            small_meta["width"] = small_meta["width"] // scale
        return smaller, small_meta


if __name__ == "__main__":
    buffer_buildings(
        path_to_esm=snakemake.input.esm,
        distance=int(snakemake.wildcards.distance),
        path_to_output=snakemake.output[0]
    )
