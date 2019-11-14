"""Rules to determine the potential of renewables on a raster map and per administrative unit.

(1) Determine technical potential on a raster map.
(2) Define administrative units.
(3) Based on scenarios, restrict the potential, and allocate it to administrative units.

"""


rule category_of_technical_eligibility:
    message:
        "Determine upper bound surface eligibility for renewables based on land cover, slope, bathymetry, and settlements."
    input:
        "src/technical_eligibility.py",
        rules.land_cover_in_europe.output,
        rules.slope_in_europe.output,
        rules.bathymetry_in_europe.output,
        rules.settlements.output.buildings,
        rules.settlements.output.urban_greens
    output:
        "build/technically-eligible-land.tif"
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT + " {CONFIG_FILE}"


rule capacityfactor_of_technical_eligibility:
    message:
        "Determine capacityfactor of eligibility category."
    input:
        "src/technically_eligible_capacityfactor.py",
        rules.category_of_technical_eligibility.output,
        expand(
            "build/capacityfactors/{technology}-time-average.tif",
            technology=["rooftop-pv", "open-field-pv", "wind-onshore", "wind-offshore"]
        )
    output:
        "build/technically-eligible-capacityfactor-pv-prio.tif",
        "build/technically-eligible-capacityfactor-wind-prio.tif"
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT + " {CONFIG_FILE}"


rule area_of_technical_eligibility:
    message:
        "Quantify the area that is technically eligible for renewables."
    input:
        "src/technically_eligible_area.py",
        rules.category_of_technical_eligibility.output,
        rules.settlements.output.buildings,
        "data/ratio-esm-available.txt"
    output:
        "build/technically-eligible-area-km2.tif"
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT


rule capacity_of_technical_eligibility:
    message:
        "Quantify the capacity that is technically eligible for renewables."
    input:
        "src/technically_eligible_capacity.py",
        rules.category_of_technical_eligibility.output,
        rules.area_of_technical_eligibility.output,
        "data/roof-statistics.csv"
    output:
        "build/technically-eligible-capacity-pv-prio-mw.tif",
        "build/technically-eligible-capacity-wind-prio-mw.tif",
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT + " {CONFIG_FILE}"


rule electricity_yield_of_technical_eligibility:
    message:
        "Quantify the max annual electricity yield that is technically eligible for renewables."
    input:
        "src/technically_eligible_electricity_yield.py",
        rules.category_of_technical_eligibility.output,
        rules.capacity_of_technical_eligibility.output,
        rules.capacityfactor_of_technical_eligibility.output
    output:
        "build/technically-eligible-electricity-yield-pv-prio-twh.tif",
        "build/technically-eligible-electricity-yield-wind-prio-twh.tif",
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT


rule buffer_buildings_high_resolution:
    message:
        "Buffer buildings of cell N{wildcards.north}E{wildcards.east} using squared distance of {wildcards.distance} m."
    input:
        src = "src/buffer_buildings.py",
        esm = "data/esm/200km_2p5m_N{north}E{east}/200km_2p5m_N{north}E{east}.TIF"
    output: temp("build/buffered-buildings-2p5m-N{north}E{east}-{distance}m.tif")
    conda: "../envs/default.yaml"
    script: "../src/buffer_buildings.py"


rule buffer_buildings_in_europe:
    message: "Merge building buffer data to study area and warp to study resolution."
    input:
        high_res = expand(
            "build/buffered-buildings-2p5m-N{north}E{east}-{{distance}}m.tif",
            north=[28, 30, 32, 34, 36],
            east=[40, 42, 44]
        ) + ["build/buffered-buildings-2p5m-N32E46-{distance}m.tif"],
        reference = rules.land_cover_in_europe.output
    output: "build/buildings-buffer-{distance}m-europe.tif"
    conda: "../envs/default.yaml"
    threads: config["snakemake"]["max-threads"]
    shadow: "minimal"
    shell:
        """
        rio merge {input.high_res} -o build/tmp-buffer.tif
        rio warp build/tmp-buffer.tif -o {output} --like {input.reference} \
        --resampling average --threads {threads}
        """


rule units:
    message: "Form units of layer {wildcards.layer} by remixing NUTS, LAU, and GADM."
    input:
        "src/units.py",
        rules.administrative_borders_nuts.output,
        rules.administrative_borders_lau.output,
        rules.administrative_borders_gadm.output
    output:
        "build/{layer}/units.geojson"
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT + " {wildcards.layer} {CONFIG_FILE}"


rule local_land_cover:
    message: "Land cover statistics per unit of layer {wildcards.layer}."
    input:
        units = rules.units.output,
        land_cover = rules.land_cover_in_europe.output,
        src = "src/geojson_to_csv.py"
    output:
        "build/{layer}/land-cover.csv"
    conda: "../envs/default.yaml"
    shell:
        """
        fio cat {input.units} | \
        rio zonalstats -r {input.land_cover} --prefix 'lc_' --categorical | \
        {PYTHON} {input.src} -a id -a lc_11 -a lc_14 -a lc_20 -a lc_30 -a lc_40 \
        -a lc_50 -a lc_60 -a lc_70 -a lc_90 -a lc_100 -a lc_110 -a lc_120 -a lc_130 \
        -a lc_140 -a lc_150 -a lc_160 -a lc_170 -a lc_180 -a lc_190 -a lc_200 \
        -a lc_210 -a lc_220 -a lc_230 > \
        {output}
        """


rule local_built_up_area:
    message: "Determine the built up area for administrative units in layer {wildcards.layer}."
    input:
        "src/built_up_area.py",
        rules.settlements.output.built_up,
        rules.units.output
    output:
        "build/{layer}/built-up-areas.csv"
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT


rule population:
    message: "Allocate population to units of layer {wildcards.layer}."
    input:
        units = rules.units.output,
        population = rules.population_in_europe.output,
        src_geojson = "src/geojson_to_csv.py",
        src_population = "src/population.py",
        land_cover = rules.local_land_cover.output
    output:
        "build/{layer}/population.csv"
    conda: "../envs/default.yaml"
    shell:
        """
        crs=$(rio info --crs {input.population})
        fio cat --dst_crs "$crs" {input.units} | \
        rio zonalstats -r {input.population} --prefix 'population_' --stats sum | \
        {PYTHON} {input.src_geojson} -a id -a population_sum -a proper | \
        {PYTHON} {input.src_population} {input.units} {input.land_cover} > \
        {output}
        """


rule demand:
    message: "Allocate electricity demand to units of layer {wildcards.layer}."
    input:
        "src/spatial_demand.py",
        rules.electricity_demand_national.output,
        rules.industry.output,
        rules.units.output,
        rules.population.output
    output:
        "build/{layer}/demand.csv"
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT


rule eez_eligibility:
    message:
        "Allocate eligible land to exclusive economic zones using {threads} threads."
    input:
        src = "src/eligibility_local.py",
        regions = rules.eez_in_europe.output,
        eligibility = rules.category_of_technical_eligibility.output
    output:
        "build/eez-eligibility.csv"
    threads: config["snakemake"]["max-threads"]
    conda: "../envs/default.yaml"
    shell:
        PYTHON + " {input.src} offshore {input.regions} {input.eligibility} {output} {threads}"


rule shared_coast:
    message: "Determine share of coast length between eez and units of layer {wildcards.layer} using {threads} threads."
    input:
        "src/shared_coast.py",
        rules.units.output,
        rules.eez_in_europe.output
    output:
        "build/{layer}/shared-coast.csv"
    threads: config["snakemake"]["max-threads"]
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT + " {threads}"


rule potentials:
    message:
        "Determine the constrained potentials for layer {wildcards.layer} in scenario {wildcards.scenario}."
    input:
        "src/potentials.py",
        rules.units.output,
        rules.eez_in_europe.output,
        rules.shared_coast.output,
        rules.electricity_yield_of_technical_eligibility.output,
        rules.category_of_technical_eligibility.output,
        rules.land_cover_in_europe.output,
        rules.protected_areas_in_europe.output
    output:
        "build/{layer}/{scenario}/potentials.csv"
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT + " {wildcards.scenario} {CONFIG_FILE}"


rule potentials_polished:
    message: "Polish potential data for publication."
    input:
        potentials = rules.potentials.output[0],
        demand = rules.demand.output[0]
    output: "build/{layer}/{scenario}/potentials-polished.csv"
    run:
        import pandas as pd

        potentials = pd.read_csv(input.potentials, index_col=0)
        demand = pd.read_csv(input.demand, index_col=0)
        potentials["Demand [TWh/yr]"] = demand["demand_twh_per_year"]
        potentials.rename(columns={
            "rooftop_pv_twh_per_year": "Roof mounted PV [TWh/yr]",
            "open_field_pv_twh_per_year": "Open field PV [TWh/yr]",
            "onshore_wind_twh_per_year": "Onshore wind [TWh/yr]",
            "offshore_wind_twh_per_year": "Offshore wind [TWh/yr]"
        }).to_csv(output[0], index=True, header=True, float_format="%.1f")


rule areas:
    message:
        "Determine eligible areas for layer {wildcards.layer} in scenario {wildcards.scenario}."
    input:
        "src/areas.py",
        rules.units.output,
        rules.eez_in_europe.output,
        rules.shared_coast.output,
        rules.area_of_technical_eligibility.output,
        rules.category_of_technical_eligibility.output,
        rules.land_cover_in_europe.output,
        rules.protected_areas_in_europe.output
    output:
        "build/{layer}/{scenario}/areas.csv"
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT + " {wildcards.scenario} {CONFIG_FILE}"


rule capacities:
    message:
        "Determine installable capacities for layer {wildcards.layer} in scenario {wildcards.scenario}."
    input:
        "src/capacities.py",
        rules.units.output,
        rules.eez_in_europe.output,
        rules.shared_coast.output,
        rules.capacity_of_technical_eligibility.output,
        rules.electricity_yield_of_technical_eligibility.output,
        rules.category_of_technical_eligibility.output,
        rules.land_cover_in_europe.output,
        rules.protected_areas_in_europe.output,
        rules.buffer_buildings_in_europe.output
    output:
        "build/{layer}/{scenario}/{distance}/capacities.csv"
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT + " {wildcards.scenario} {CONFIG_FILE}"


rule normed_potentials:
    message:
        "Determine potentials relative to demand for layer {wildcards.layer} "
        "in scenario {wildcards.scenario}."
    input:
        "src/potentials_normed.py",
        rules.demand.output,
        rules.potentials.output
    output:
        "build/{layer}/{scenario}/normed-potentials.csv"
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT


rule footprint:
    message: "Determine the land footprint of the renewable potential for layer {wildcards.layer} "
             "in scenario {wildcards.scenario}."
    input:
        "src/footprint.py",
        rules.category_of_technical_eligibility.output,
        rules.area_of_technical_eligibility.output,
        rules.electricity_yield_of_technical_eligibility.output,
        rules.land_cover_in_europe.output,
        rules.protected_areas_in_europe.output,
        rules.units.output
    output:
        "build/{layer}/{scenario}/footprint.csv"
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT + " {wildcards.scenario} {CONFIG_FILE}"


rule necessary_land:
    message: "Determine the necessary potentials to become autarkic of layer {wildcards.layer} "
             "in scenario {wildcards.scenario} given rooftop PV share {wildcards.pvshare}%."
    input:
        "src/necessary_land.py",
        rules.demand.output,
        rules.potentials.output,
        rules.footprint.output,
        rules.local_built_up_area.output,
        rules.units.output
    output:
        "build/{layer}/{scenario}/necessary-land-when-pv-{pvshare}%.csv"
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT + " {wildcards.pvshare}"


rule scenario_results:
    message: "Merge all results of scenario {wildcards.scenario} for layer {wildcards.layer}."
    input:
        units = rules.units.output,
        demand = rules.demand.output,
        population = rules.population.output,
        constrained_potentials = rules.potentials.output,
        normed_potentials = rules.normed_potentials.output
    output:
        "build/{layer}/{scenario}/merged-results.gpkg"
    run:
        import pandas as pd
        import geopandas as gpd

        gpd.read_file(input.units[0]).merge(
            pd.concat(
                [pd.read_csv(path, index_col="id") for path in input[1:]],
                axis=1
            ),
            left_on="id",
            right_index=True
        ).to_file(output[0], driver="GPKG")
