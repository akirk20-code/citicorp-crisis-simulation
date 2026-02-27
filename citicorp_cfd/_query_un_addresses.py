#!/usr/bin/env python3
"""Query building addresses in the UN HQ area from NYC Open Data."""

import json
import sys

try:
    import requests
except ImportError:
    print("Installing requests...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "requests", "-q"])
    import requests

# UN HQ center
CENTER_LAT = 40.7489
CENTER_LON = -73.9680
RADIUS = 350  # meters, slightly larger to capture edges

M_PER_DEG_LAT = 111024.0
M_PER_DEG_LON = 84352.0
FT_TO_M = 0.3048

# Bounding box
nw_lat = CENTER_LAT + RADIUS / M_PER_DEG_LAT
nw_lon = CENTER_LON - RADIUS / M_PER_DEG_LON
se_lat = CENTER_LAT - RADIUS / M_PER_DEG_LAT
se_lon = CENTER_LON + RADIUS / M_PER_DEG_LON

# ----- Step 1: Get building footprints from Socrata GeoJSON -----
print("Querying NYC Building Footprints (Socrata)...")
where = f"within_box(the_geom, {nw_lat}, {nw_lon}, {se_lat}, {se_lon})"
params = {"$where": where, "$limit": 500}
resp = requests.get(
    "https://data.cityofnewyork.us/resource/5zhs-2jue.geojson",
    params=params, timeout=60
)
resp.raise_for_status()
data = resp.json()
features = data.get("features", [])
print(f"  Received {len(features)} features")

# Filter to >= 5m and sort by height
buildings = []
for f in features:
    props = f.get("properties", {})
    h_ft = float(props.get("height_roof") or 0)
    h_m = h_ft * FT_TO_M
    if h_m >= 5.0:
        buildings.append({
            "bin": str(props.get("bin") or ""),
            "base_bbl": str(props.get("base_bbl") or ""),
            "mappluto_bbl": str(props.get("mappluto_bbl") or ""),
            "name": props.get("name") or "",
            "height_m": h_m,
            "year": str(props.get("construction_year") or ""),
            "feat_code": str(props.get("feature_code") or ""),
            "ground_elev_ft": float(props.get("ground_elevation") or 0),
        })

buildings.sort(key=lambda b: b["height_m"], reverse=True)
print(f"  Buildings >= 5m: {len(buildings)}")

# ----- Step 2: Get addresses from PLUTO via mappluto_bbl -----
print("\nQuerying MapPLUTO for addresses...")
# Use mappluto_bbl (matches PLUTO's bbl field, which is numeric)
bbls = list(set(b["mappluto_bbl"] for b in buildings if b["mappluto_bbl"]))
print(f"  Unique mappluto_bbls to look up: {len(bbls)}")

addresses = {}  # keyed by mappluto_bbl string

# Query in chunks using numeric comparison (PLUTO bbl is numeric)
for i in range(0, len(bbls), 15):
    chunk = bbls[i:i+15]
    bbl_filter = " OR ".join([f"bbl={b}" for b in chunk])
    try:
        pluto_resp = requests.get(
            "https://data.cityofnewyork.us/resource/64uk-42ks.json",
            params={
                "$where": bbl_filter,
                "$select": "bbl,address,zipcode,ownername,numfloors,yearbuilt,landuse,bldgclass",
                "$limit": 200
            },
            timeout=30
        )
        if pluto_resp.status_code == 200:
            for p in pluto_resp.json():
                # PLUTO bbl is like "1013357501.00000000" — normalize to int string
                raw_bbl = p.get("bbl", "")
                try:
                    norm_bbl = str(int(float(raw_bbl)))
                except (ValueError, TypeError):
                    norm_bbl = raw_bbl
                addresses[norm_bbl] = p
        else:
            print(f"  PLUTO chunk {i}: HTTP {pluto_resp.status_code}")
    except Exception as e:
        print(f"  PLUTO chunk {i} error: {e}")

print(f"  Got PLUTO data for {len(addresses)} BBLs")

# ----- Step 3: Print results -----
print(f"\n{'='*110}")
print(f"UN HQ Area Building Inventory ({len(buildings)} buildings >= 5m)")
print(f"Center: {CENTER_LAT}, {CENTER_LON} | Radius: {RADIUS}m")
print(f"{'='*110}")
print(f"{'#':>3} {'BIN':>7} {'H(m)':>6} {'Yr':>5} {'Fl':>3} {'Address':<50} {'Zip':<6} {'Owner':<30}")
print(f"{'-'*110}")

for i, b in enumerate(buildings):
    mbb = b["mappluto_bbl"]
    bin_val = b["bin"]

    pluto = addresses.get(mbb, {})
    addr = pluto.get("address", "(no address)")
    zipcode = pluto.get("zipcode", "")
    floors_raw = pluto.get("numfloors", "")
    owner = pluto.get("ownername", "")
    try:
        floors = str(int(float(floors_raw))) if floors_raw else ""
    except (ValueError, TypeError):
        floors = ""

    name = b["name"]
    if name:
        addr = f"{addr} [{name}]"

    # Truncate owner
    if len(owner) > 28:
        owner = owner[:26] + ".."

    print(f"{i+1:3d} {bin_val:>7} {b['height_m']:>6.1f} {b['year']:>5} {floors:>3} {addr:<50} {zipcode:<6} {owner:<30}")

print(f"\n{'='*110}")
