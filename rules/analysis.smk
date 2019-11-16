"""Rules analysing the scenario results."""


wildcard_constraints:
    plot_suffix = "((png)|(tif))" # can plot tif or png


rule necessary_land_overview:
    message: "Create table showing the fraction of land needed to become autarkic in scenario {wildcards.scenario} "
             "for rooftop PV share {wildcards.pvshare}%."
    input:
        nec_land = expand("build/{layer}/{{scenario}}/necessary-land-when-pv-{{pvshare}}%.csv", layer=config["layers"].keys()),
        demand = expand("build/{layer}/demand.csv", layer=config["layers"].keys())
    output:
        "build/{scenario}/overview-necessary-land-when-pv-{pvshare}%.csv"
    run:
        import pandas as pd
        nec_lands = [pd.read_csv(path, index_col=0)["fraction_non_built_up_land_necessary"] for path in input.nec_land]
        nec_roofs = [pd.read_csv(path, index_col=0)["fraction_roofs_necessary"] for path in input.nec_land]
        roof_pv_gens = [pd.read_csv(path, index_col=0)["rooftop_pv_generation_twh_per_year"] for path in input.nec_land]
        demands = [pd.read_csv(path, index_col=0)["demand_twh_per_year"] for path in input.demand]
        roof_pv_shares = [roof_pv_gen / demand for roof_pv_gen, demand in zip(roof_pv_gens, demands)]
        data = pd.DataFrame(
            index=[name.capitalize() for name in config["layers"].keys()],
            data={
                "Average land use [%]": [nec_land.mean() * 100 for nec_land in nec_lands],
                "Average roof space use [%]": [nec_roof.mean() * 100 for nec_roof in nec_roofs],
                "Average roof-mounted PV share [%]:": [roof_pv_share.mean() * 100 for roof_pv_share in roof_pv_shares]
            }
        )
        data.index.name = "Level"
        data.to_csv(
            output[0],
            index=True,
            header=True,
            float_format="%.0f"
        )


rule normed_potential_boxplots:
    message: "Plot ranges of relative potential for scenario {wildcards.scenario}."
    input:
        "src/vis/potentials_normed_boxplot.py",
        "build/municipal/{scenario}/merged-results.gpkg"
    output:
        "build/{scenario}/normed-potentials-boxplots.{plot_suffix}"
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT


rule potentials_sufficiency_map:
    message: "Plot potential sufficiency maps for scenario {wildcards.scenario}."
    input:
        "src/vis/potentials_sufficiency_map.py",
        "build/continental/{scenario}/merged-results.gpkg",
        "build/national/{scenario}/merged-results.gpkg",
        "build/regional/{scenario}/merged-results.gpkg",
        "build/municipal/{scenario}/merged-results.gpkg"
    output:
        "build/{scenario}/sufficient-potentials-map.{plot_suffix}"
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT


rule exclusion_layers_plot:
    message: "Visualise the exclusion layers for {wildcards.country_code}."
    input:
        "src/vis/exclusion_layers.py",
        "build/national/units.geojson",
        rules.land_cover_in_europe.output,
        rules.slope_in_europe.output,
        rules.protected_areas_in_europe.output,
        rules.settlements.output.buildings,
    output:
        "build/exclusion-layers-{country_code}.{plot_suffix}"
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT + " {wildcards.country_code}"


rule layer_overview:
    message: "Overview table of administrative layers."
    input:
        expand("build/{layer}/units.geojson", layer=config["layers"].keys())
    output:
        "build/overview-administrative-levels.csv"
    run:
        import pandas as pd
        import geopandas as gpd

        def format_source(source_name):
            source = source_name[:-1].upper()
            if source in ["NUTS", "LAU"]:
                source = source + " [@eurostat:2015]"
            elif source == "GADM":
                source = source + " [@GADM:2018]"
            return source

        layer_names = [path_to_file.split("/")[1] for path_to_file in input]
        sources = [[format_source(source) for source in set(config["layers"][layer_name].values())]
                   for layer_name in layer_names]
        sources = [", ".join(set(sources_per_layer)) for sources_per_layer in sources]
        number_units = [len(gpd.read_file(path_to_file).index) for path_to_file in input]

        pd.DataFrame({
            "Level": [name.capitalize() for name in layer_names],
            "Number units": number_units,
            "Source of shape data": sources
        }).to_csv(output[0], index=False, header=True)


rule scenario_overview:
    message: "Overview over scenario on all levels."
    input:
        normed_potentials = expand("build/{layer}/{{scenario}}/normed-potentials.csv", layer=config["layers"].keys()),
        population = expand("build/{layer}/population.csv", layer=config["layers"].keys())
    output:
        "build/{scenario}/overview.csv"
    run:
        import pandas as pd

        def add_layer(df, path_to_file):
            layer_name = path_to_file.split("/")[1]
            df["layer"] = layer_name
            return df

        ALL_LAYERS = ["continental", "national", "regional", "municipal"]

        normed_potential = pd.concat([pd.read_csv(path, index_col="id").pipe(lambda df: add_layer(df, path))
                                     for path in input.normed_potentials])
        population = pd.concat([pd.read_csv(path, index_col="id").pipe(lambda df: add_layer(df, path))
                               for path in input.population])

        affected = normed_potential.normed_potential < 1.0
        high_density = population.density_p_per_km2 > 1000
        overview = pd.DataFrame(
            index=[path_to_file.split("/")[1] for path_to_file in input.normed_potentials]
        )
        overview["units affected [%]"] = (normed_potential.loc[affected].groupby("layer").count() /
                                          normed_potential.groupby("layer").count() * 100).reindex(ALL_LAYERS).iloc[:, 0]
        overview["of which dense units [%]"] = (normed_potential.loc[affected & high_density].groupby("layer").count() /
                                                normed_potential.loc[affected].groupby("layer").count() * 100).reindex(ALL_LAYERS).iloc[:, 0]
        overview["people affected [%]"] = (population.loc[affected].groupby("layer").population_sum.sum() /
                                           population.groupby("layer").population_sum.sum() * 100).reindex(ALL_LAYERS)
        overview["of which from dense units [%]"] = (population.loc[affected & high_density].groupby("layer").population_sum.sum() /
                                                     population.loc[affected].groupby("layer").population_sum.sum() * 100).reindex(ALL_LAYERS)

        overview.fillna(0).transpose().to_csv(output[0], index=True, header=True, float_format="%.1f")


rule wind_capacity_per_distance_plot:
    message: "Visualise the dependency between wind capacity and distance."
    input:
        "src/vis/wind_capacity.py",
        results = expand(
            "build/{layer}/{scenario}/{distance}/capacities.csv",
            scenario=["technical-potential", "technical-potential-env-protection"],
            layer=["national"],
            distance=[600, 800, 1000, 1200]
        )
    output:
        "build/wind-capacity-per-distance.png"
    conda: "../envs/default.yaml"
    script: "../src/vis/wind_capacity.py"


rule wind_capacity_per_distance_map:
    message: "Visualise the dependency between wind capacity and distance on a map."
    input:
        "src/vis/wind_capacity_map.py",
        results = expand(
            "build/{layer}/{scenario}/{distance}/capacities.csv",
            scenario=["technical-potential-env-protection"],
            layer=["regional", "municipal"],
            distance=[600, 1000]
        ),
        regions = "build/regional/units.geojson",
        municipalities = "build/municipal/units.geojson"
    output:
        germany = "build/wind-capacity-per-distance-map-germany.png",
        brandenburg = "build/wind-capacity-per-distance-map-brandenburg.png"
    conda: "../envs/default.yaml"
    script: "../src/vis/wind_capacity_map.py"
