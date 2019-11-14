import math
import gc

import numpy as np
import rasterio
import scipy.ndimage


BUILDING_CODE = 50
RESOLUTION_IN_M = 2.5
SCALING_FACTOR = 10 # results in data 100 times smaller


def buffer_buildings(path_to_esm, distance, path_to_output):
    with rasterio.open(path_to_esm) as src:
        esm = src.read(1)
        meta = src.meta

    # TODO should use circles, but scipy function is extremely memory hungry
    buffered = buffer_buildings_using_squares(esm, buffer_in_m=distance)
    del esm
    gc.collect() # immediately remove unnecessary data to avoid memory peaks

    smaller, small_meta = downsample(buffered, meta, SCALING_FACTOR)
    del buffered
    gc.collect()

    with rasterio.open(path_to_output, 'w', **small_meta) as dst:
        dst.write(smaller)


def buffer_buildings_using_circles(esm, buffer_in_m):
    mask = esm == BUILDING_CODE
    filter_mask = circle_mask(math.ceil(buffer_in_m / RESOLUTION_IN_M))
    # split array, because operation is too memory hungry
    left, right = tuple(np.hsplit(mask, 2))
    upper_left, lower_left = tuple(np.vsplit(left, 2))
    upper_right, lower_right = tuple(np.vsplit(right, 2))
    del left, right
    # compute buffers
    upper_left = buffer_parts_of_buildings_using_circles(upper_left, filter_mask)
    upper_right = buffer_parts_of_buildings_using_circles(upper_right, filter_mask)
    lower_left = buffer_parts_of_buildings_using_circles(lower_left, filter_mask)
    lower_right = buffer_parts_of_buildings_using_circles(lower_right, filter_mask)
    # restack subarrays
    left = np.vstack((upper_left, lower_left))
    right = np.vstack((upper_right, lower_right))
    return np.hstack((left, right))


def buffer_parts_of_buildings_using_circles(esm, filter_mask):
    return scipy.ndimage.maximum_filter(
        esm,
        footprint=filter_mask,
        mode='constant',
        cval=0
    )


def buffer_buildings_using_squares(esm, buffer_in_m):
    mask = esm == BUILDING_CODE
    return scipy.ndimage.maximum_filter(
        mask,
        size=math.ceil(buffer_in_m / RESOLUTION_IN_M) * 2 + 1,
        mode='constant',
        cval=0
    )


# adapted from https://stackoverflow.com/a/18354475/1856079
def circle_mask(radius):
    """
    Return a boolean mask for a circle.
    """
    matrix_size = radius * 2 + 1
    x, y = np.ogrid[:matrix_size, :matrix_size]
    cx, cy = radius, radius

    # convert cartesian --> polar coordinates
    r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy)

    # circular mask
    circmask = r2 <= radius * radius

    return circmask * 1


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
