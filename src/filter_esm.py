import gc

import numpy as np
import rasterio
import scipy.ndimage


def filter_esm(path_to_esm_tile, path_to_output):
    with rasterio.open(path_to_esm_tile, "r") as src:
        esm = src.read(1)
        meta = src.meta

    filtered = scipy.ndimage.median_filter(
        esm,
        size=7,
        mode='constant',
        cval=0
    )
    del esm
    gc.collect() # immediately remove unnecessary data to avoid memory peaks

    meta["compress"] = "lzw"
    with rasterio.open(path_to_output, 'w', **meta) as dst:
        dst.write(filtered, 1)


def building_filter(window):
    # Functionwise, this is a great filter for the problems I see.
    # Unfortunately, it's roughly 40 times slower than the median filter.
    center = window[12] # ASSUME window size 5
    uniques = np.unique(window)
    if center == 50:
        if uniques[0] == 50: # only buildings here
            return center
        elif uniques[-2] >= 30:
            return center
        else:
            return 20 # anything but building, doesn't matter
    else:
        return center


if __name__ == "__main__":
    filter_esm(
        path_to_esm_tile=snakemake.input.esm,
        path_to_output=snakemake.output[0]
    )
