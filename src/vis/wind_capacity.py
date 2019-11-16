import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

RED = "#A01914"
BLUE = "#4F6DB8"


def plot_wind_capacity_per_distance(paths_to_results, path_to_plot):
    sns.set_context('paper')
    data = pd.DataFrame.from_records(
        columns=["potential_type", "distance_m", "onshore_potential"],
        data=[
            (pottype(path_to_result), distance(path_to_result), potential(path_to_result))
            for path_to_result in paths_to_results
        ]
    )
    data["onshore_potential"] = data.groupby("potential_type").onshore_potential.transform(lambda x: x / x.max())

    fig = plt.figure(figsize=(8, 4), constrained_layout=True)
    ax = fig.subplots(1, 1)
    sns.barplot(
        data=data,
        x="distance_m",
        y="onshore_potential",
        hue="potential_type",
        palette=[RED, BLUE],
        ax=ax
    )
    fig.savefig(path_to_plot, dpi=600, transparent=False)


def pottype(path_to_result):
    potential = path_to_result.split("/")[-3]
    if potential == "technical-potential":
        return "ignoring environmental protection"
    else:
        return "considering environmental protection"


def distance(path_to_result):
    return int(path_to_result.split("/")[-2])


def potential(path_to_result):
    return pd.read_csv(path_to_result).loc[0, "onshore_wind_mw"]


if __name__ == "__main__":
    plot_wind_capacity_per_distance(
        paths_to_results=snakemake.input.results,
        path_to_plot=snakemake.output[0]
    )
