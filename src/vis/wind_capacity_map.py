import pandas as pd
import geopandas as gpd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt


RED = "#A01914"
BLUE = "#4F6DB8"
GREY = "#C0C0C0"
RED_TO_BLUE = [ # from https://gka.github.io using lightness correction
    '#002d6e', '#375aa2', '#6f8ad1', '#a7bffa',
    '#f5f5f5', '#fdad97', '#e36b55', '#b23125', '#720000'
]
CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("signature-BlRd", RED_TO_BLUE)
NORM = matplotlib.colors.Normalize(vmin=0, vmax=1)
PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"
EDGE_WIDTH = 0.06
EDGE_COLOR = "white"

EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "

MAP_MIN_X = 4050000
MAP_MIN_Y = 2650000
MAP_MAX_X = 4700000
MAP_MAX_Y = 3650000


def plot_wind_capacity_per_distance_map(paths_to_results, path_to_regional_shapes,
                                        path_to_municipal_shapes, path_to_plot):
    sns.set_context('paper')
    regions, municipalities = read_results(paths_to_results, path_to_regional_shapes, path_to_municipal_shapes)

    fig = plt.figure(figsize=(8, 5), constrained_layout=True)
    axes = fig.subplots(1, 2).flatten()
    plot_layer(regions, "Bundesländer", "a", 0.06, axes[0])
    plot_layer(municipalities, "Gemeinden", "b", 0.0125, axes[1])
    plot_colorbar(fig, axes)
    fig.savefig(path_to_plot, dpi=600, transparent=False)


def read_results(paths_to_results, path_to_region_shapes, path_to_municipal_shapes):
    regional_results = filter(lambda x: "regional" in x, paths_to_results)
    municipal_results = filter(lambda x: "municipal" in x, paths_to_results)
    regions = gpd.read_file(path_to_region_shapes).to_crs(EPSG_3035_PROJ4).set_index("id")
    regions = regions.assign(
        **{distance(path_to_result): pd.read_csv(path_to_result).set_index("id").loc[:, "onshore_wind_mw"]
           for path_to_result in regional_results}
    )
    regions["reduced_potential"] = 1 - regions["1000"] / regions["600"]

    municipalities = gpd.read_file(path_to_municipal_shapes).to_crs(EPSG_3035_PROJ4).set_index("id")
    municipalities = municipalities.assign(
        **{distance(path_to_result): pd.read_csv(path_to_result).set_index("id").loc[:, "onshore_wind_mw"]
           for path_to_result in municipal_results}
    )
    municipalities["reduced_potential"] = 1 - municipalities["1000"] / municipalities["600"]
    return regions, municipalities


def pottype(path_to_result):
    potential = path_to_result.split("/")[-3]
    if potential == "technical-potential":
        return "ignoring environmental protection"
    else:
        return "considering environmental protection"


def distance(path_to_result):
    return path_to_result.split("/")[-2]


def potential(path_to_result):
    pot = pd.read_csv(path_to_result).loc[:, "onshore_wind_mw"]
    if pottype(path_to_result) == "considering environmental protection":
        return pot * 10
    else:
        return pot


def plot_layer(units, layer_name, panel_id, edge_width, ax):
    ax.set_aspect('equal')
    invalids = units[units.reduced_potential.isna()]
    units[~units.isin(invalids)].plot(
        linewidth=edge_width,
        edgecolor=EDGE_COLOR,
        column="reduced_potential",
        vmin=NORM.vmin,
        vmax=NORM.vmax,
        cmap=CMAP,
        ax=ax
    )
    if not invalids.empty:
        invalids.plot(
            linewidth=edge_width,
            edgecolor=EDGE_COLOR,
            facecolor=GREY,
            ax=ax
        )
    ax.set_xlim(MAP_MIN_X, MAP_MAX_X)
    ax.set_ylim(MAP_MIN_Y, MAP_MAX_Y)
    ax.set_xticks([])
    ax.set_yticks([])
    sns.despine(ax=ax, top=True, bottom=True, left=True, right=True)
    ax.annotate(layer_name, xy=[0.2, 0.95], xycoords='axes fraction')
    ax.annotate(panel_id, xy=[0.05, 0.95], xycoords='axes fraction',
                fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)


def plot_colorbar(fig, axes):
    s_m = matplotlib.cm.ScalarMappable(cmap=CMAP, norm=NORM)
    cbar = fig.colorbar(s_m, ax=axes, fraction=1, aspect=35, shrink=0.65)
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1])
    cbar.set_ticklabels(["0%", "25%", "50%", "75%", "100%"])
    cbar.outline.set_linewidth(0)
    cbar.ax.set_ylabel('Reduktion des Windpotenzials bei \n1000m Abstand im Vergleich zu 600m', rotation=90)


if __name__ == "__main__":
    plot_wind_capacity_per_distance_map(
        paths_to_results=snakemake.input.results,
        path_to_regional_shapes=snakemake.input.regions,
        path_to_municipal_shapes=snakemake.input.municipalities,
        path_to_plot=snakemake.output[0]
    )
