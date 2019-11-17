import pandas as pd
import geopandas as gpd
import shapely
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt


RED = "#A01914"
BLUE = "#4F6DB8"
GREY = "#C0C0C0"
BLUE_TO_RED = [ # from https://gka.github.io using lightness correction
    '#002d6e', '#375aa2', '#6f8ad1', '#a7bffa',
    '#f5f5f5', '#fdad97', '#e36b55', '#b23125', '#720000'
]
BLUE_TO_RED.reverse()
CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("signature-BlRd", BLUE_TO_RED)
NORM = matplotlib.colors.Normalize(vmin=0, vmax=1)
PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"
EDGE_WIDTH = 0.06
EDGE_COLOR = "white"

EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "

GERMANY_BOUNDING_BOX = {
    "min_x": 4050000,
    "min_y": 2650000,
    "max_x": 4700000,
    "max_y": 3650000,
}

BRANDENBURG_BOUNDING_BOX = {
    "min_x": 4400000,
    "min_y": 3100000,
    "max_x": 4700000,
    "max_y": 3450000,
}

BRANDENBURG_ID = "DEU.4_1"


def plot_wind_capacity_per_distance_map(paths_to_results, path_to_regional_shapes, path_to_municipal_shapes,
                                        path_to_germany_plot, path_to_brandenburg_plot):
    sns.set_context('paper')
    regions, municipalities = read_results(paths_to_results, path_to_regional_shapes, path_to_municipal_shapes)
    brandenburg, brandenburg_municipalities = filter_brandenburg(regions, municipalities)

    fig = plot_germany(regions, municipalities)
    fig.savefig(path_to_germany_plot, dpi=600, transparent=False)
    fig = plot_brandenburg(brandenburg, brandenburg_municipalities)
    fig.savefig(path_to_brandenburg_plot, dpi=600, transparent=False)


def plot_germany(regions, municipalities):
    fig = plt.figure(figsize=(8, 5), constrained_layout=True)
    axes = fig.subplots(1, 2).flatten()
    plot_layer(regions, "a - Bundesländer", 0.06, GERMANY_BOUNDING_BOX, axes[0])
    plot_layer(municipalities, "b - Gemeinden", 0.0125, GERMANY_BOUNDING_BOX, axes[1])
    plot_colorbar(fig, axes)
    return fig


def plot_brandenburg(regions, municipalities):
    fig = plt.figure(figsize=(8, 4), constrained_layout=True)
    axes = fig.subplots(1, 2).flatten()
    plot_layer(regions, "a – Brandenburg", 0.06, BRANDENBURG_BOUNDING_BOX, axes[0])
    plot_layer(municipalities, "b - Gemeinden", 0.0125, BRANDENBURG_BOUNDING_BOX, axes[1])
    plot_colorbar(fig, axes)
    return fig


def read_results(paths_to_results, path_to_region_shapes, path_to_municipal_shapes):
    regional_results = filter(lambda x: "regional" in x, paths_to_results)
    municipal_results = filter(lambda x: "municipal" in x, paths_to_results)
    regions = gpd.read_file(path_to_region_shapes).to_crs(EPSG_3035_PROJ4).set_index("id")
    regions = regions.assign(
        **{distance(path_to_result): pd.read_csv(path_to_result).set_index("id").loc[:, "onshore_wind_mw"]
           for path_to_result in regional_results}
    )
    regions["reduced_potential"] = regions["1000"] / regions["600"]

    municipalities = gpd.read_file(path_to_municipal_shapes).to_crs(EPSG_3035_PROJ4).set_index("id")
    municipalities = municipalities.assign(
        **{distance(path_to_result): pd.read_csv(path_to_result).set_index("id").loc[:, "onshore_wind_mw"]
           for path_to_result in municipal_results}
    )
    municipalities["reduced_potential"] = municipalities["1000"] / municipalities["600"]
    return regions, municipalities


def filter_brandenburg(regions, municipalities):
    brandenburg = regions.loc[BRANDENBURG_ID].geometry.buffer(1000)
    brandenburg.prep = shapely.prepared.prep(brandenburg)
    mask = (brandenburg.contains(gem) for gem in municipalities.geometry)
    return regions[regions.index == BRANDENBURG_ID], municipalities.loc[mask]


def pottype(path_to_result):
    potential = path_to_result.split("/")[-3]
    if potential == "technical-potential":
        return "ignoring environmental protection"
    else:
        return "considering environmental protection"


def distance(path_to_result):
    return path_to_result.split("/")[-2]


def plot_layer(units, panel_id, edge_width, bounding_box, ax):
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
    ax.set_xlim(bounding_box["min_x"], bounding_box["max_x"])
    ax.set_ylim(bounding_box["min_y"], bounding_box["max_y"])
    ax.set_xticks([])
    ax.set_yticks([])
    sns.despine(ax=ax, top=True, bottom=True, left=True, right=True)
    ax.annotate(panel_id, xy=[0.05, 0.95], xycoords='axes fraction',
                fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)


def plot_colorbar(fig, axes):
    s_m = matplotlib.cm.ScalarMappable(cmap=CMAP, norm=NORM)
    cbar = fig.colorbar(s_m, ax=axes, fraction=1, aspect=35, shrink=0.65)
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1])
    cbar.set_ticklabels(["0%", "25%", "50%", "75%", "100%"])
    cbar.outline.set_linewidth(0)
    cbar.ax.set_ylabel('Verbleibendes Windpotenzials bei \n1000m Abstand im Vergleich zu 600m', rotation=90)


if __name__ == "__main__":
    plot_wind_capacity_per_distance_map(
        paths_to_results=snakemake.input.results,
        path_to_regional_shapes=snakemake.input.regions,
        path_to_municipal_shapes=snakemake.input.municipalities,
        path_to_germany_plot=snakemake.output.germany,
        path_to_brandenburg_plot=snakemake.output.brandenburg
    )
