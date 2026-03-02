#!/usr/bin/env python3
"""Check if NYC Building Footprints API has 1978 historical data."""

import requests
import json

API_URL = "https://data.cityofnewyork.us/resource/5zhs-2jue.geojson"

# Query for Citicorp area
CITICORP_LAT = 40.7579
CITICORP_LON = -73.9690
LAT_MIN = 40.7507
LAT_MAX = 40.7651
LON_MIN = -73.9762
LON_MAX = -73.9618

print("Querying NYC API for available fields...")

# Get just a few buildings to see what fields are available
params = {
    "$where": f"within_box(the_geom, {LAT_MAX}, {LON_MIN}, {LAT_MIN}, {LON_MAX})",
    "$limit": 5,
}

try:
    resp = requests.get(API_URL, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    if data.get("features"):
        sample = data["features"][0]
        props = sample.get("properties", {})

        print("\nAvailable fields in NYC Building Footprints API:")
        print("=" * 70)
        for key, value in sorted(props.items()):
            print(f"  {key:25s}: {value}")

        # Check for construction year fields
        year_fields = [k for k in props.keys() if 'year' in k.lower() or 'cnstrct' in k.lower() or 'built' in k.lower() or 'date' in k.lower()]

        print("\n" + "=" * 70)
        if year_fields:
            print("YEAR/DATE FIELDS FOUND:")
            for field in year_fields:
                print(f"  ✓ {field}: {props.get(field)}")
        else:
            print("✗ NO CONSTRUCTION YEAR FIELDS FOUND")
            print("\nNYC Building Footprints API does NOT include construction dates.")
            print("This API only has current geometry and physical attributes.")

        print("\n" + "=" * 70)
        print("\nAlternative data sources for 1978 historical skyline:")
        print("  1. NYC PLUTO dataset (Property Land Use Tax)")
        print("     - URL: https://data.cityofnewyork.us/City-Government/Primary-Land-Use-Tax-Lot-Output-PLUTO/64uk-42ks")
        print("     - Has 'YearBuilt' field")
        print("     - Can filter: YearBuilt <= 1978")
        print("     - Caveat: Tax lot polygons, not building footprints")
        print()
        print("  2. NYC 3D Building Model (Historical)")
        print("     - URL: https://www1.nyc.gov/site/doitt/initiatives/3d-building.page")
        print("     - May have construction year metadata")
        print()
        print("  3. Emporis or CoStar databases")
        print("     - Commercial building databases with construction dates")
        print("     - Requires subscription")
        print()
        print("  4. Historical aerial photos + manual digitization")
        print("     - NYC Municipal Archives aerial photos from 1970s")
        print("     - URL: https://nycma.lunaimaging.com/luna/servlet")
        print("     - Labor-intensive but most accurate")

    else:
        print("No buildings returned from API")

except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
