"""Determine potential of renewable electricity in each administrative unit.

* Take the (only technically restricted) raster data potentials,
* add restrictions based on scenario definitions,
* allocate the onshore potentials to the administrative units,
* allocate the offshore potentials to exclusive economic zones (EEZ),
* allocate the offshore potential of EEZ to units based on the fraction of shared coast.

This is in analogy to `potentials.py` and `areas.py` but for installable capacities [MW]
rather than areas [km2] or [TWh/a].
"""
import click
import numpy as np
import pandas as pd
import rasterio
import fiona

from src.technical_eligibility import Eligibility, FOREST, FARM, OTHER
from src.potentials import Potential, decide_between_pv_and_wind, potentials_per_shape, ProtectedArea
from src.utils import Config


@click.command()
@click.argument("path_to_units")
@click.argument("path_to_eez")
@click.argument("path_to_shared_coast")
@click.argument("path_to_capacities_pv_prio")
@click.argument("path_to_capacities_wind_prio")
@click.argument("path_to_electricity_yield_pv_prio")
@click.argument("path_to_electricity_yield_wind_prio")
@click.argument("path_to_eligibility_categories")
@click.argument("path_to_land_cover")
@click.argument("path_to_protected_areas")
@click.argument("path_to_building_buffer")
@click.argument("path_to_result")
@click.argument("scenario")
@click.argument("config", type=Config())
def potentials(path_to_units, path_to_eez, path_to_shared_coast,
               path_to_capacities_pv_prio, path_to_capacities_wind_prio,
               path_to_electricity_yield_pv_prio, path_to_electricity_yield_wind_prio,
               path_to_eligibility_categories, path_to_land_cover, path_to_protected_areas,
               path_to_building_buffer, path_to_result, scenario, config):
    """Determine potential of renewable electricity in each administrative unit.

    * Take the (only technically restricted) raster data potentials,
    * add restrictions based on scenario definitions,
    * allocate the onshore potentials to the administrative units,
    * allocate the offshore potentials to exclusive economic zones (EEZ),
    * allocate the offshore potential of EEZ to units based on the fraction of shared coast.
    """
    with rasterio.open(path_to_eligibility_categories, "r") as src:
        eligibility_categories = src.read(1)
    with rasterio.open(path_to_capacities_pv_prio, "r") as src:
        transform = src.transform
        capacities_pv_prio = src.read(1)
    with rasterio.open(path_to_capacities_wind_prio, "r") as src:
        capacities_wind_prio = src.read(1)
    with rasterio.open(path_to_electricity_yield_pv_prio, "r") as src:
        transform = src.transform
        electricity_yield_pv_prio = src.read(1)
    with rasterio.open(path_to_electricity_yield_wind_prio, "r") as src:
        electricity_yield_wind_prio = src.read(1)
    with rasterio.open(path_to_land_cover, "r") as src:
        land_cover = src.read(1)
    with rasterio.open(path_to_protected_areas, "r") as src:
        protected_areas = src.read(1)
    with rasterio.open(path_to_building_buffer, "r") as src:
        building_buffer = 1 - src.read(1)
    with fiona.open(path_to_units, "r") as src:
        unit_ids = [feature["properties"]["id"] for feature in src]
        unit_geometries = [feature["geometry"] for feature in src]
    with fiona.open(path_to_eez, "r") as src:
        eez_ids = [feature["properties"]["id"] for feature in src]
        eez_geometries = [feature["geometry"] for feature in src]
    shared_coasts = pd.read_csv(path_to_shared_coast, index_col=0)

    capacities_pv_prio, capacities_wind_prio = apply_scenario_config(
        potential_pv_prio=capacities_pv_prio,
        potential_wind_prio=capacities_wind_prio,
        categories=eligibility_categories,
        land_cover=land_cover,
        protected_areas=protected_areas,
        building_buffer=building_buffer,
        scenario_config=config["scenarios"][scenario]
    )
    electricity_yield_pv_prio, electricity_yield_wind_prio = apply_scenario_config(
        potential_pv_prio=electricity_yield_pv_prio,
        potential_wind_prio=electricity_yield_wind_prio,
        categories=eligibility_categories,
        land_cover=land_cover,
        protected_areas=protected_areas,
        building_buffer=building_buffer,
        scenario_config=config["scenarios"][scenario]
    )
    capacities_pv_prio, capacities_wind_prio = decide_between_pv_and_wind(
        potential_pv_prio=capacities_pv_prio,
        potential_wind_prio=capacities_wind_prio,
        electricity_yield_pv_prio=electricity_yield_pv_prio,
        electricity_yield_wind_prio=electricity_yield_wind_prio,
        eligibility_categories=eligibility_categories
    )

    onshore_potentials = pd.DataFrame(
        index=unit_ids,
        data={
            potential.capacity_name: potentials_per_shape(
                eligibilities=potential.eligible_on,
                potential_map=(capacities_pv_prio if "pv" in str(potential).lower()
                               else capacities_wind_prio),
                eligibility_categories=eligibility_categories,
                shapes=unit_geometries,
                transform=transform
            )
            for potential in Potential.onshore()
        }
    )
    offshore_eez_potentials = pd.DataFrame(
        index=eez_ids,
        data={
            potential.capacity_name: potentials_per_shape(
                eligibilities=potential.eligible_on,
                potential_map=(capacities_pv_prio if "pv" in str(potential).lower()
                               else capacities_wind_prio),
                eligibility_categories=eligibility_categories,
                shapes=eez_geometries,
                transform=transform
            )
            for potential in Potential.offshore()
        }
    )
    offshore_potentials = pd.DataFrame(
        data=shared_coasts.dot(offshore_eez_potentials),
        columns=[potential.capacity_name for potential in Potential.offshore()]
    )
    potentials = pd.concat([onshore_potentials, offshore_potentials], axis=1)
    potentials.index.name = "id"
    potentials.to_csv(
        path_to_result,
        header=True,
        index=True
    )


def apply_scenario_config(potential_pv_prio, potential_wind_prio, categories,
                          land_cover, protected_areas, building_buffer, scenario_config):
    """Limit potential in each pixel based on scenario config."""

    # share-rooftops-used
    share_rooftops_used = scenario_config["share-rooftops-used"]
    mask = categories == Eligibility.ROOFTOP_PV
    potential_pv_prio[mask] = potential_pv_prio[mask] * share_rooftops_used
    potential_wind_prio[mask] = potential_wind_prio[mask] * share_rooftops_used

    # share-forest-used-for-wind
    share_forest_used_for_wind = scenario_config["share-forest-used-for-wind"]
    mask = np.isin(land_cover, FOREST) & (categories != Eligibility.ROOFTOP_PV)
    potential_pv_prio[mask] = potential_pv_prio[mask] * share_forest_used_for_wind
    potential_wind_prio[mask] = potential_wind_prio[mask] * share_forest_used_for_wind

    # share-other-land-used
    share_other_land_used = scenario_config["share-other-land-used"]
    mask = np.isin(land_cover, OTHER) & (categories != Eligibility.ROOFTOP_PV)
    potential_pv_prio[mask] = potential_pv_prio[mask] * share_other_land_used
    potential_wind_prio[mask] = potential_wind_prio[mask] * share_other_land_used

    # share-farmland-used
    share_farmland_used = scenario_config["share-farmland-used"]
    mask = np.isin(land_cover, FARM) & (categories != Eligibility.ROOFTOP_PV)
    potential_pv_prio[mask] = potential_pv_prio[mask] * share_farmland_used
    potential_wind_prio[mask] = potential_wind_prio[mask] * share_farmland_used

    # share-offshore-used
    share_offshore_used = scenario_config["share-offshore-used"]
    mask = categories == Eligibility.OFFSHORE_WIND
    potential_pv_prio[mask] = potential_pv_prio[mask] * share_offshore_used
    potential_wind_prio[mask] = potential_wind_prio[mask] * share_offshore_used

    # pv-on-farmland
    pv_on_farmland = scenario_config["pv-on-farmland"]
    if not pv_on_farmland:
        mask = np.isin(land_cover, FARM) & (categories == Eligibility.ONSHORE_WIND_AND_PV)
        potential_pv_prio[mask] = 0

    # share-protected-areas-used
    use_protected_areas = scenario_config["use-protected-areas"]
    if not use_protected_areas:
        mask = (protected_areas == ProtectedArea.PROTECTED) & (categories != Eligibility.ROOFTOP_PV)
        potential_pv_prio[mask] = 0
        potential_wind_prio[mask] = 0

    # don't build wind close to buildings
    mask = (categories == Eligibility.ONSHORE_WIND)
    potential_pv_prio[mask] = potential_pv_prio[mask] * building_buffer[mask]
    potential_wind_prio[mask] = potential_wind_prio[mask] * building_buffer[mask]
    mask = (categories == Eligibility.ONSHORE_WIND_AND_PV)
    potential_wind_prio[mask] = potential_wind_prio[mask] * building_buffer[mask]

    return potential_pv_prio, potential_wind_prio


if __name__ == "__main__":
    potentials()
