#!/usr/bin/env python3
"""Analyze building heights around Citicorp Center from NYC Open Data API."""

import requests
import statistics

# Citicorp Center location
CITICORP_LAT = 40.7579
CITICORP_LON = -73.9690

# Query parameters (same as generate_stl.py)
API_URL = "https://data.cityofnewyork.us/resource/5zhs-2jue.geojson"
FT_TO_M = 0.3048

# Bounding box (slightly larger than CFD domain for context)
# Domain: ±360m = ±0.003° lat, ±0.005° lon
LAT_MIN = 40.7507  # -360m
LAT_MAX = 40.7651  # +360m
LON_MIN = -73.9762 # -360m
LON_MAX = -73.9618 # +360m

print("Querying NYC Building Footprints API...")
print(f"Bounding box: {LAT_MIN:.4f} to {LAT_MAX:.4f}°N, {LON_MIN:.4f} to {LON_MAX:.4f}°W")

# SoQL query: within_box(the_geom, top_lat, left_lon, bottom_lat, right_lon)
params = {
    "$where": f"within_box(the_geom, {LAT_MAX}, {LON_MIN}, {LAT_MIN}, {LON_MAX})",
    "$limit": 10000,
    "$select": "bin,name,height_roof,ground_elevation"
}

try:
    resp = requests.get(API_URL, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    buildings = data.get("features", [])
    print(f"\nReceived {len(buildings)} buildings from API\n")

    # Extract heights
    heights_ft = []
    heights_m = []
    citicorp_bin = None
    tallest_neighbors = []

    for feature in buildings:
        props = feature.get("properties", {})
        bin_str = props.get("bin", "")
        name = props.get("name", "Unknown")
        height_ft = props.get("height_roof")

        if height_ft is not None:
            height_ft = float(height_ft)
            height_m = height_ft * FT_TO_M

            # Check if this is Citicorp
            if bin_str in ("1035879", "1087931"):
                citicorp_bin = bin_str
                print(f"✓ Found Citicorp Center: BIN {bin_str}")
                print(f"  Reported height: {height_ft:.1f} ft = {height_m:.1f} m")
                print(f"  (Note: API data may not include crown/roof detail)\n")
                continue

            # Collect non-Citicorp buildings
            heights_ft.append(height_ft)
            heights_m.append(height_m)

            # Track tallest neighbors
            if height_m > 100:  # Buildings over 100m
                tallest_neighbors.append((name or f"BIN {bin_str}", height_m, bin_str))

    # Sort neighbors by height
    tallest_neighbors.sort(key=lambda x: x[1], reverse=True)

    # Statistics
    if heights_m:
        print("=" * 70)
        print(f"BUILDING HEIGHT STATISTICS (Excluding Citicorp)")
        print("=" * 70)
        print(f"Total buildings analyzed: {len(heights_m)}")
        print(f"\nHeight range:")
        print(f"  Minimum:  {min(heights_m):6.1f} m ({min(heights_ft):5.0f} ft)")
        print(f"  Maximum:  {max(heights_m):6.1f} m ({max(heights_ft):5.0f} ft)")
        print(f"  Mean:     {statistics.mean(heights_m):6.1f} m ({statistics.mean(heights_ft):5.0f} ft)")
        print(f"  Median:   {statistics.median(heights_m):6.1f} m ({statistics.median(heights_ft):5.0f} ft)")
        print(f"  Std dev:  {statistics.stdev(heights_m):6.1f} m")

        # Height distribution
        print(f"\nHeight distribution:")
        ranges = [
            (0, 20, "Very low (0-20m)"),
            (20, 50, "Low-rise (20-50m)"),
            (50, 100, "Mid-rise (50-100m)"),
            (100, 150, "High-rise (100-150m)"),
            (150, 200, "Tall (150-200m)"),
            (200, 250, "Very tall (200-250m)"),
            (250, 300, "Supertall (250-300m)"),
        ]

        for low, high, label in ranges:
            count = sum(1 for h in heights_m if low <= h < high)
            pct = count / len(heights_m) * 100
            bar = "█" * int(pct / 2)
            print(f"  {label:25s}: {count:4d} buildings ({pct:5.1f}%) {bar}")

        # Tallest neighbors
        print(f"\nTOP 20 TALLEST NEIGHBORS:")
        print(f"{'Rank':<6} {'Height':>12} {'BIN':>10}  {'Name':<40}")
        print("-" * 70)
        for i, (name, height_m, bin_str) in enumerate(tallest_neighbors[:20], 1):
            height_ft = height_m / FT_TO_M
            print(f"{i:<6} {height_m:6.1f}m ({height_ft:4.0f}ft) {bin_str:>10}  {name[:40]}")

        # Comparison to Citicorp
        citicorp_height = 278.9  # meters (actual with crown)
        print(f"\n" + "=" * 70)
        print(f"CITICORP COMPARISON")
        print("=" * 70)
        print(f"Citicorp height (with crown): {citicorp_height:.1f} m (915 ft)")
        print(f"Tallest neighbor:             {max(heights_m):.1f} m ({max(heights_ft):.0f} ft)")
        print(f"Citicorp / Tallest ratio:     {citicorp_height / max(heights_m):.2f}×")
        print(f"\nNeighbors taller than 200m:   {sum(1 for h in heights_m if h > 200)}")
        print(f"Neighbors taller than 150m:   {sum(1 for h in heights_m if h > 150)}")
        print(f"Neighbors taller than 100m:   {sum(1 for h in heights_m if h > 100)}")

        # Height percentile
        sorted_heights = sorted(heights_m)
        citicorp_percentile = sum(1 for h in sorted_heights if h < citicorp_height) / len(sorted_heights) * 100
        print(f"\nCiticorp is taller than {citicorp_percentile:.1f}% of surrounding buildings")

    else:
        print("No building height data found!")

except Exception as e:
    print(f"Error querying API: {e}")
    import traceback
    traceback.print_exc()
