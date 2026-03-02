#!/usr/bin/env python3
"""Query NYC buildings constructed before 1978."""

import requests

API_URL = "https://data.cityofnewyork.us/resource/5zhs-2jue.geojson"
FT_TO_M = 0.3048

LAT_MIN, LAT_MAX = 40.7507, 40.7651
LON_MIN, LON_MAX = -73.9762, -73.9618

print("Checking if construction_year field is populated...")

# Query for buildings WITH construction year data
params = {
    "$where": f"within_box(the_geom, {LAT_MAX}, {LON_MIN}, {LAT_MIN}, {LON_MAX}) AND construction_year IS NOT NULL",
    "$limit": 100,
    "$select": "bin,name,height_roof,construction_year",
}

try:
    resp = requests.get(API_URL, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    buildings = data.get("features", [])
    print(f"\nBuildings WITH construction_year data: {len(buildings)}")

    if buildings:
        print("\nSample buildings with construction year:")
        for i, feature in enumerate(buildings[:10], 1):
            props = feature.get("properties", {})
            name = props.get("name", "Unknown")
            year = props.get("construction_year")
            height_ft = props.get("height_roof", 0)
            height_m = float(height_ft) * FT_TO_M if height_ft else 0
            print(f"  {i}. {name[:40]:40s} | {year} | {height_m:6.1f}m")

        # Try filtering to pre-1978
        print("\nAttempting to filter to pre-1978 buildings...")
        params_1978 = {
            "$where": f"within_box(the_geom, {LAT_MAX}, {LON_MIN}, {LAT_MIN}, {LON_MAX}) AND construction_year <= 1978",
            "$limit": 10000,
        }
        resp2 = requests.get(API_URL, params=params_1978, timeout=30)
        resp2.raise_for_status()
        data2 = resp2.json()
        pre_1978 = data2.get("features", [])
        print(f"Buildings with construction_year <= 1978: {len(pre_1978)}")

        if pre_1978:
            print("\nYES - 1978 historical data IS AVAILABLE!")
            print(f"Found {len(pre_1978)} buildings from 1978 or earlier")

            # Statistics
            years = [f.get("properties", {}).get("construction_year") for f in pre_1978]
            years = [y for y in years if y]
            if years:
                print(f"\nConstruction year range: {min(years)} - {max(years)}")

    else:
        print("\nNO - construction_year field is NOT populated in this dataset")
        print("\nAlternative: Use NYC PLUTO dataset")
        print("URL: https://data.cityofnewyork.us/City-Government/Primary-Land-Use-Tax-Lot-Output-PLUTO/64uk-42ks")
        print("Field: YearBuilt")

except Exception as e:
    print(f"Error: {e}")
