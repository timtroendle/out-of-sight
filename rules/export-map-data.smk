"""Rules to generate data that is used for the online map."""

WEB_MERCATOR = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 "\
               "+units=m +nadgrids=@null +wktext  +no_defs"

localrules: export_map_data


rule export_map_data:
    message: "Export everything."
    input:
        "build/export/map/national-boundaries-national-level.shp",
        "build/export/map/national-boundaries-regional-level.shp",
        "build/export/map/national-boundaries-municipal-level.shp",
        "build/export/map/national--technical-potential-env-protection-all-distances.mbtiles",
        "build/export/map/regional--technical-potential-env-protection-all-distances.mbtiles",
        "build/export/map/municipal--technical-potential-env-protection-all-distances.mbtiles",


rule national_boundaries:
    message: "Export national boundaries on {wildcards.layer} level."
    input: "build/{layer}/units.geojson"
    output: "build/export/map/national-boundaries-{layer}-level.shp"
    run:
        import geopandas as gpd

        units = gpd.read_file(input[0])[["country_code", "geometry"]]
        units.dissolve(by="country_code").geometry.to_crs(WEB_MERCATOR).to_file(output[0])


rule scenario_result:
    message: "Export data of {wildcards.layer} layer and all distances."
    input:
        potential_600 = "build/{layer}/technical-potential-env-protection/600/capacities.csv",
        potential_800 = "build/{layer}/technical-potential-env-protection/800/capacities.csv",
        potential_1000 = "build/{layer}/technical-potential-env-protection/1000/capacities.csv",
        potential_1200 = "build/{layer}/technical-potential-env-protection/1200/capacities.csv",
        units = "build/{layer}/units.geojson"
    output: "build/export/map/{layer}--technical-potential-env-protection-all-distances.geojson"
    run:
        import pandas as pd
        import geopandas as gpd

        def read_wind_potential(path_to_potential):
            pot = (
                pd
                .read_csv(path_to_potential)
                .rename(columns={"id": "unit_id"})
                .set_index("unit_id")
                .loc[:, "onshore_wind_mw"]
            )
            pot[pot < 1] = 0 # neglect potentials below 1 MW
            return pot / 8 * 10 # from 8 MW/km2 to 10 MW/km2

        potential_1200 = read_wind_potential(input.potential_1200)
        potential_1000 = read_wind_potential(input.potential_1000)
        potential_800 = read_wind_potential(input.potential_800)
        potential_600 = read_wind_potential(input.potential_600)

        potential = pd.DataFrame(
            index=potential_600.index,
            data={
                "onshore_wind_mw_600": potential_600,
                "onshore_wind_mw_800": potential_800,
                "onshore_wind_mw_1000": potential_1000,
                "onshore_wind_mw_1200": potential_1200,
                "onshore_wind_rel_1000_to_600": potential_1000 / potential_600
            }
        )

        (
            gpd
            .read_file(input.units)
            .merge(potential, left_on="id", right_index=True)
            .replace("Germany", "Deutschland")
            .to_file(output[0], driver="GeoJSON")
        )


def zoom_range_parameters(wildcards):
    filename = wildcards.filename
    if "continental" in filename:
        return "-z4"
    elif "national" in filename:
        return "-Z3 -z6"
    elif "regional" in filename:
        return "-Z5 -z9"
    elif "municipal" in filename:
        return "-Z9 -z12"
    else:
        return "-zg" # determine automatically


rule tiles:
    message: "Create tiles from '{wildcards.filename}.geojson'."
    input: "build/export/map/{filename}.geojson"
    params: zoom_range = zoom_range_parameters
    output: "build/export/map/{filename}.mbtiles"
    conda: "../envs/webmap.yaml"
    shell: "tippecanoe -o {output} {params.zoom_range} --generate-ids {input}"
