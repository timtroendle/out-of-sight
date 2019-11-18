import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

EXPANSION_MW_2019 = 287
BLUE = "#4F6DB8"
COLOR = BLUE


def capacity_expansion(path_to_data, path_to_plot):
    wind_deployment = pd.read_csv(path_to_data)
    wind_deployment = pd.concat([
        wind_deployment,
        pd.DataFrame(index=["Jahre", "Zubau", "Kumuliert"], data=["2019*", EXPANSION_MW_2019, pd.np.nan]).T
    ])

    sns.set_context('paper')
    fig = plt.figure(figsize=(8, 3), constrained_layout=True)
    ax = fig.subplots(1, 1)

    sns.barplot(
        data=wind_deployment,
        x="Jahre",
        y="Zubau",
        palette=[COLOR],
        ax=ax
    )
    ax.set_xlabel("Jahr")
    ax.set_ylabel("Zubau (MW)")
    sns.despine()
    fig.savefig(path_to_plot, dpi=600, transparent=False)


if __name__ == "__main__":
    capacity_expansion(
        path_to_data=snakemake.input.winddata,
        path_to_plot=snakemake.output[0]
    )
