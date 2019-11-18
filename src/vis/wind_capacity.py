import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

RED = "#A01914"
BLUE = "#4F6DB8"


def plot_wind_capacity_per_distance(paths_to_results, path_to_plot):
    sns.set_context('paper')
    data = pd.DataFrame.from_records(
        columns=["potential_type", "distance_m", "onshore_potential_gw"],
        data=[
            (pottype(path_to_result), distance(path_to_result), potential_gw(path_to_result))
            for path_to_result in paths_to_results
        ]
    )
    data2 = data.copy()
    data3 = data.copy()
    data2["onshore_potential_gw"] = data2["onshore_potential_gw"] / 8 * 10
    data3["onshore_potential_gw"] = data3["onshore_potential_gw"] / 8 * 12

    fig = plt.figure(figsize=(8, 4), constrained_layout=True)
    ax = fig.subplots(1, 1)
    sns.barplot(
        data=pd.concat([data, data2, data3]),
        x="distance_m",
        y="onshore_potential_gw",
        palette=[BLUE],
        errwidth=4,
        ax=ax
    )
    sns.despine()
    ax.set_xlabel("Mindestabstand (m)")
    ax.set_ylabel("Technisches Potenzial (GW)")
    fig.savefig(path_to_plot, dpi=600, transparent=False)


def pottype(path_to_result):
    potential = path_to_result.split("/")[-3]
    if potential == "technical-potential":
        return "ignoring environmental protection"
    else:
        return "considering environmental protection"


def distance(path_to_result):
    return int(path_to_result.split("/")[-2])


def potential_gw(path_to_result):
    return pd.read_csv(path_to_result).loc[0, "onshore_wind_mw"] / 1000


if __name__ == "__main__":
    plot_wind_capacity_per_distance(
        paths_to_results=snakemake.input.results,
        path_to_plot=snakemake.output[0]
    )
