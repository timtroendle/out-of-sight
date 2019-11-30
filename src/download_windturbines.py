import requests
import json


def download(path_to_file):
    overpass_url = "http://overpass-api.de/api/interpreter"
    overpass_query = """
    [out:json];
    area["ISO3166-1"="DE"][admin_level=2];
    node["generator:source"="wind"]["power"="generator"](area);
    out center;
    """
    response = requests.get(overpass_url,
                            params={'data': overpass_query})
    turbines = response.json()

    with open(path_to_file, "w") as dst:
        json.dump(turbines, dst)


if __name__ == "__main__":
    download(snakemake.output[0])
