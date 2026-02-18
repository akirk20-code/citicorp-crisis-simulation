#!/usr/bin/env python3
"""
generate_stl_osm.py

Queries OpenStreetMap via Overpass API to extract building footprints and heights
around Citicorp Center (40.7579, -73.9690), then generates STL files for CFD simulation.

Outputs:
  - citicorp_tower.stl: Main tower with 45° slanted roof and 45° rotation
  - citicorp_stilts.stl: 10-story stilts at corners
  - surroundings_osm.stl: All other buildings from OSM data

Dependencies: requests (can be installed via: pip install requests)
Fallback: Uses only standard library if requests unavailable
"""

import sys
import json
import math
from typing import List, Tuple, Dict, Optional

# Try to import requests, fallback to urllib if unavailable
try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    import urllib.request
    import urllib.parse
    HAS_REQUESTS = False
    print("Warning: 'requests' library not found, using urllib fallback")


# =============================================================================
# CONFIGURATION
# =============================================================================

# Citicorp Center coordinates (WGS84)
CITICORP_LAT = 40.7579
CITICORP_LON = -73.9690

# Domain size: ±360m from center
DOMAIN_RADIUS = 360.0  # meters

# OSM Overpass API endpoint
OVERPASS_URL = "https://overpass-api.de/api/interpreter"

# Height estimation parameters
DEFAULT_FLOOR_HEIGHT = 3.5  # meters per floor
CITICORP_HEIGHT = 279.0  # meters (59 floors × 4.73m/floor)
CITICORP_ROOF_ANGLE = 45.0  # degrees
CITICORP_ROTATION = 45.0  # degrees from north
CITICORP_STILT_STORIES = 10
CITICORP_STILT_HEIGHT = CITICORP_STILT_STORIES * DEFAULT_FLOOR_HEIGHT

# Coordinate conversion (approximate for NYC latitude)
METERS_PER_DEGREE_LAT = 111320.0
METERS_PER_DEGREE_LON = 111320.0 * math.cos(math.radians(CITICORP_LAT))


# =============================================================================
# COORDINATE CONVERSION
# =============================================================================

def latlon_to_meters(lat: float, lon: float) -> Tuple[float, float]:
    """
    Convert WGS84 lat/lon to local Cartesian coordinates (meters).
    Origin at Citicorp Center.

    Returns: (x, y) where x is East, y is North
    """
    x = (lon - CITICORP_LON) * METERS_PER_DEGREE_LON
    y = (lat - CITICORP_LAT) * METERS_PER_DEGREE_LAT
    return x, y


def compute_bounding_box() -> Tuple[float, float, float, float]:
    """
    Compute bounding box in WGS84 coordinates.

    Returns: (south, west, north, east)
    """
    delta_lat = DOMAIN_RADIUS / METERS_PER_DEGREE_LAT
    delta_lon = DOMAIN_RADIUS / METERS_PER_DEGREE_LON

    south = CITICORP_LAT - delta_lat
    north = CITICORP_LAT + delta_lat
    west = CITICORP_LON - delta_lon
    east = CITICORP_LON + delta_lon

    return south, west, north, east


# =============================================================================
# OVERPASS API QUERY
# =============================================================================

def query_overpass(bbox: Tuple[float, float, float, float]) -> Optional[Dict]:
    """
    Query Overpass API for buildings in bounding box.

    Args:
        bbox: (south, west, north, east) in WGS84 degrees

    Returns:
        JSON response as dict, or None on error
    """
    south, west, north, east = bbox

    # Overpass QL query
    # Note: bbox format is [south,west,north,east]
    query = f"""
    [bbox:{south},{west},{north},{east}];
    (
      way["building"];
      relation["building"];
    );
    out geom;
    """

    print(f"Querying Overpass API for bbox: {bbox}")
    print(f"  South: {south:.6f}, West: {west:.6f}")
    print(f"  North: {north:.6f}, East: {east:.6f}")

    try:
        if HAS_REQUESTS:
            response = requests.post(OVERPASS_URL, data={"data": query}, timeout=60)
            response.raise_for_status()
            return response.json()
        else:
            # Fallback using urllib
            data = urllib.parse.urlencode({"data": query}).encode("utf-8")
            req = urllib.request.Request(OVERPASS_URL, data=data)
            with urllib.request.urlopen(req, timeout=60) as response:
                return json.loads(response.read().decode("utf-8"))

    except Exception as e:
        print(f"Error querying Overpass API: {e}")
        return None


# =============================================================================
# HEIGHT ESTIMATION FROM OSM TAGS
# =============================================================================

def estimate_height(tags: Dict[str, str]) -> float:
    """
    Estimate building height from OSM tags.

    Priority:
      1. height tag (meters)
      2. building:levels × 3.5m
      3. building:height
      4. roof:height
      5. Default: 15m (approx 4-5 stories)

    Args:
        tags: OSM tags dictionary

    Returns:
        Estimated height in meters
    """
    # Try explicit height tag
    if "height" in tags:
        try:
            h = float(tags["height"].replace("m", "").strip())
            if h > 0:
                return h
        except ValueError:
            pass

    # Try building:height
    if "building:height" in tags:
        try:
            h = float(tags["building:height"].replace("m", "").strip())
            if h > 0:
                return h
        except ValueError:
            pass

    # Try building:levels
    if "building:levels" in tags:
        try:
            levels = float(tags["building:levels"])
            if levels > 0:
                return levels * DEFAULT_FLOOR_HEIGHT
        except ValueError:
            pass

    # Try roof:height as fallback
    if "roof:height" in tags:
        try:
            h = float(tags["roof:height"].replace("m", "").strip())
            if h > 0:
                return h
        except ValueError:
            pass

    # Default: 4-5 story building
    return 15.0


def is_citicorp_building(footprint: List[Tuple[float, float]], tags: Dict[str, str]) -> bool:
    """
    Determine if a building is Citicorp Center.

    Heuristics:
      - Centroid near (0, 0) in local coordinates
      - Height > 200m or has "Citicorp" in name
      - Large footprint area

    Args:
        footprint: List of (x, y) coordinates in meters
        tags: OSM tags

    Returns:
        True if likely Citicorp Center
    """
    # Check name tags
    name = tags.get("name", "").lower()
    if "citicorp" in name or "601 lexington" in name:
        return True

    # Check centroid proximity (within 50m of origin)
    centroid_x = sum(x for x, y in footprint) / len(footprint)
    centroid_y = sum(y for x, y in footprint) / len(footprint)
    dist = math.sqrt(centroid_x**2 + centroid_y**2)

    if dist < 50.0:
        # Check height
        height = estimate_height(tags)
        if height > 200.0:
            return True

        # Check footprint area (Citicorp is roughly 65m × 65m)
        area = compute_polygon_area(footprint)
        if area > 3000.0:  # ~3000-4500 m²
            return True

    return False


def compute_polygon_area(footprint: List[Tuple[float, float]]) -> float:
    """
    Compute polygon area using shoelace formula.

    Args:
        footprint: List of (x, y) coordinates

    Returns:
        Area in square meters
    """
    if len(footprint) < 3:
        return 0.0

    area = 0.0
    n = len(footprint)
    for i in range(n):
        j = (i + 1) % n
        area += footprint[i][0] * footprint[j][1]
        area -= footprint[j][0] * footprint[i][1]

    return abs(area) / 2.0


# =============================================================================
# STL GEOMETRY GENERATION
# =============================================================================

class STLWriter:
    """
    Simple ASCII STL file writer.
    """

    def __init__(self, filename: str):
        self.filename = filename
        self.facets = []

    def add_triangle(self, v1: Tuple[float, float, float],
                     v2: Tuple[float, float, float],
                     v3: Tuple[float, float, float]):
        """
        Add a triangle to the STL.
        Vertices should be in counter-clockwise order when viewed from outside.

        Args:
            v1, v2, v3: (x, y, z) coordinates
        """
        # Compute normal vector (right-hand rule)
        ax, ay, az = v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]
        bx, by, bz = v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2]

        nx = ay * bz - az * by
        ny = az * bx - ax * bz
        nz = ax * by - ay * bx

        # Normalize
        length = math.sqrt(nx**2 + ny**2 + nz**2)
        if length > 0:
            nx /= length
            ny /= length
            nz /= length

        self.facets.append((nx, ny, nz, v1, v2, v3))

    def write(self):
        """Write STL file to disk."""
        print(f"Writing {len(self.facets)} facets to {self.filename}")

        with open(self.filename, "w") as f:
            f.write(f"solid geometry\n")

            for nx, ny, nz, v1, v2, v3 in self.facets:
                f.write(f"  facet normal {nx:.6e} {ny:.6e} {nz:.6e}\n")
                f.write(f"    outer loop\n")
                f.write(f"      vertex {v1[0]:.6e} {v1[1]:.6e} {v1[2]:.6e}\n")
                f.write(f"      vertex {v2[0]:.6e} {v2[1]:.6e} {v2[2]:.6e}\n")
                f.write(f"      vertex {v3[0]:.6e} {v3[1]:.6e} {v3[2]:.6e}\n")
                f.write(f"    endloop\n")
                f.write(f"  endfacet\n")

            f.write(f"endsolid geometry\n")


def extrude_footprint(footprint: List[Tuple[float, float]], height: float) -> List[Tuple]:
    """
    Extrude a 2D footprint to 3D with flat roof.

    Args:
        footprint: List of (x, y) in meters
        height: Building height in meters

    Returns:
        List of triangles as (v1, v2, v3) tuples
    """
    triangles = []
    n = len(footprint)

    if n < 3:
        return triangles

    # Bottom face (z=0) - triangulate using fan method
    for i in range(1, n - 1):
        v1 = (footprint[0][0], footprint[0][1], 0.0)
        v2 = (footprint[i][0], footprint[i][1], 0.0)
        v3 = (footprint[i+1][0], footprint[i+1][1], 0.0)
        triangles.append((v1, v3, v2))  # Reversed for downward normal

    # Top face (z=height)
    for i in range(1, n - 1):
        v1 = (footprint[0][0], footprint[0][1], height)
        v2 = (footprint[i][0], footprint[i][1], height)
        v3 = (footprint[i+1][0], footprint[i+1][1], height)
        triangles.append((v1, v2, v3))

    # Side walls
    for i in range(n):
        j = (i + 1) % n

        x1, y1 = footprint[i]
        x2, y2 = footprint[j]

        # Two triangles per wall segment
        v1 = (x1, y1, 0.0)
        v2 = (x2, y2, 0.0)
        v3 = (x2, y2, height)
        v4 = (x1, y1, height)

        triangles.append((v1, v2, v3))
        triangles.append((v1, v3, v4))

    return triangles


def generate_citicorp_tower() -> STLWriter:
    """
    Generate hand-crafted Citicorp tower with 45° slanted roof and 45° rotation.

    Returns:
        STLWriter with tower geometry
    """
    stl = STLWriter("citicorp_tower.stl")

    # Tower dimensions (approximate)
    width = 60.0  # meters (roughly 65m in reality)
    height_base = CITICORP_HEIGHT - 10.0  # Base of slanted roof
    height_peak = CITICORP_HEIGHT  # Peak of roof

    # Rotation: 45° from north (counterclockwise about z-axis)
    theta = math.radians(CITICORP_ROTATION)
    cos_t = math.cos(theta)
    sin_t = math.sin(theta)

    def rotate(x: float, y: float) -> Tuple[float, float]:
        """Rotate point about origin."""
        return x * cos_t - y * sin_t, x * sin_t + y * cos_t

    # Base footprint (square centered at origin, then rotated)
    half_w = width / 2.0
    corners = [
        rotate(-half_w, -half_w),
        rotate(half_w, -half_w),
        rotate(half_w, half_w),
        rotate(-half_w, half_w),
    ]

    # Extrude main tower body (flat top at height_base)
    for tri in extrude_footprint(corners, height_base):
        stl.add_triangle(*tri)

    # Slanted roof: 45° slope from south to north
    # Peak along north edge, base along south edge
    # After rotation, this creates the characteristic slant

    # Roof coordinates (before rotation)
    roof_south_low = [
        (-half_w, -half_w, height_base),
        (half_w, -half_w, height_base),
    ]
    roof_north_high = [
        (-half_w, half_w, height_peak),
        (half_w, half_w, height_peak),
    ]

    # Rotate roof corners
    roof_corners_low = [
        (rotate(x, y)[0], rotate(x, y)[1], z) for x, y, z in roof_south_low
    ]
    roof_corners_high = [
        (rotate(x, y)[0], rotate(x, y)[1], z) for x, y, z in roof_north_high
    ]

    # Roof faces (4 triangles forming slanted plane)
    stl.add_triangle(roof_corners_low[0], roof_corners_low[1], roof_corners_high[1])
    stl.add_triangle(roof_corners_low[0], roof_corners_high[1], roof_corners_high[0])

    # Roof side walls (triangular gables)
    # West gable
    stl.add_triangle(
        roof_corners_low[0],
        roof_corners_high[0],
        (roof_corners_low[0][0], roof_corners_low[0][1], height_base)
    )
    # East gable
    stl.add_triangle(
        roof_corners_low[1],
        (roof_corners_low[1][0], roof_corners_low[1][1], height_base),
        roof_corners_high[1]
    )
    # South wall (vertical)
    stl.add_triangle(
        corners[0] + (height_base,),
        corners[1] + (height_base,),
        roof_corners_low[1]
    )
    stl.add_triangle(
        corners[0] + (height_base,),
        roof_corners_low[1],
        roof_corners_low[0]
    )
    # North wall (vertical at peak)
    stl.add_triangle(
        corners[2] + (height_base,),
        roof_corners_high[0],
        roof_corners_high[1]
    )
    stl.add_triangle(
        corners[2] + (height_base,),
        roof_corners_high[1],
        corners[3] + (height_base,)
    )

    return stl


def generate_citicorp_stilts() -> STLWriter:
    """
    Generate Citicorp stilts (10-story columns at corners).

    Returns:
        STLWriter with stilt geometry
    """
    stl = STLWriter("citicorp_stilts.stl")

    # Stilt dimensions
    stilt_width = 8.0  # meters (approximate column width)
    tower_width = 60.0

    # Tower rotation
    theta = math.radians(CITICORP_ROTATION)
    cos_t = math.cos(theta)
    sin_t = math.sin(theta)

    def rotate(x: float, y: float) -> Tuple[float, float]:
        return x * cos_t - y * sin_t, x * sin_t + y * cos_t

    # Four corner positions (before rotation)
    half_tw = tower_width / 2.0
    corner_offsets = [
        (-half_tw, -half_tw),
        (half_tw, -half_tw),
        (half_tw, half_tw),
        (-half_tw, half_tw),
    ]

    # Generate stilts at each corner
    for offset_x, offset_y in corner_offsets:
        # Stilt footprint (small square)
        half_sw = stilt_width / 2.0
        stilt_corners = [
            rotate(offset_x - half_sw, offset_y - half_sw),
            rotate(offset_x + half_sw, offset_y - half_sw),
            rotate(offset_x + half_sw, offset_y + half_sw),
            rotate(offset_x - half_sw, offset_y + half_sw),
        ]

        # Extrude stilt from ground to stilt height
        for tri in extrude_footprint(stilt_corners, CITICORP_STILT_HEIGHT):
            stl.add_triangle(*tri)

    return stl


def generate_surrounding_buildings(osm_data: Dict) -> STLWriter:
    """
    Generate STL for all surrounding buildings from OSM data.

    Args:
        osm_data: Overpass API JSON response

    Returns:
        STLWriter with surrounding building geometry
    """
    stl = STLWriter("surroundings_osm.stl")

    if not osm_data or "elements" not in osm_data:
        print("Warning: No OSM data available for surrounding buildings")
        return stl

    building_count = 0
    skipped_count = 0

    for element in osm_data["elements"]:
        # Skip non-building elements
        if element.get("type") not in ["way", "relation"]:
            continue

        tags = element.get("tags", {})
        if "building" not in tags:
            continue

        # Extract footprint coordinates
        footprint_latlon = []

        if element["type"] == "way":
            # Way: nodes contain geometry
            if "geometry" in element:
                footprint_latlon = [(node["lat"], node["lon"]) for node in element["geometry"]]
            else:
                skipped_count += 1
                continue

        elif element["type"] == "relation":
            # Relation: extract outer way (for multipolygons)
            members = element.get("members", [])
            for member in members:
                if member.get("role") == "outer" and "geometry" in member:
                    footprint_latlon = [(node["lat"], node["lon"]) for node in member["geometry"]]
                    break

            if not footprint_latlon:
                skipped_count += 1
                continue

        # Convert to local Cartesian coordinates
        footprint = [latlon_to_meters(lat, lon) for lat, lon in footprint_latlon]

        # Skip if outside domain
        in_domain = False
        for x, y in footprint:
            if abs(x) <= DOMAIN_RADIUS and abs(y) <= DOMAIN_RADIUS:
                in_domain = True
                break

        if not in_domain:
            skipped_count += 1
            continue

        # Skip Citicorp itself
        if is_citicorp_building(footprint, tags):
            print(f"Skipping Citicorp building: {tags.get('name', 'unnamed')}")
            skipped_count += 1
            continue

        # Estimate height
        height = estimate_height(tags)

        # Generate building geometry
        for tri in extrude_footprint(footprint, height):
            stl.add_triangle(*tri)

        building_count += 1

    print(f"Generated {building_count} surrounding buildings")
    print(f"Skipped {skipped_count} elements")

    return stl


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """Main execution function."""
    print("="*70)
    print("Citicorp Center STL Generator from OpenStreetMap")
    print("="*70)
    print()

    # Step 1: Compute bounding box
    bbox = compute_bounding_box()
    print(f"Domain: {DOMAIN_RADIUS}m radius around Citicorp Center")
    print(f"Center: ({CITICORP_LAT:.6f}, {CITICORP_LON:.6f})")
    print()

    # Step 2: Query OSM data
    print("Querying OpenStreetMap Overpass API...")
    osm_data = query_overpass(bbox)

    if osm_data:
        num_elements = len(osm_data.get("elements", []))
        print(f"Retrieved {num_elements} elements from OSM")
    else:
        print("Warning: Failed to retrieve OSM data")
        print("Continuing with hand-crafted Citicorp geometry only")

    print()

    # Step 3: Generate Citicorp tower
    print("Generating Citicorp tower (hand-crafted geometry)...")
    tower_stl = generate_citicorp_tower()
    tower_stl.write()
    print()

    # Step 4: Generate Citicorp stilts
    print("Generating Citicorp stilts...")
    stilts_stl = generate_citicorp_stilts()
    stilts_stl.write()
    print()

    # Step 5: Generate surrounding buildings
    print("Generating surrounding buildings from OSM...")
    surroundings_stl = generate_surrounding_buildings(osm_data)
    surroundings_stl.write()
    print()

    # Summary
    print("="*70)
    print("STL Generation Complete")
    print("="*70)
    print(f"  citicorp_tower.stl:       {len(tower_stl.facets)} facets")
    print(f"  citicorp_stilts.stl:      {len(stilts_stl.facets)} facets")
    print(f"  surroundings_osm.stl:     {len(surroundings_stl.facets)} facets")
    print()
    print("Next steps:")
    print("  1. Verify STL files with viewer (e.g., MeshLab, Blender)")
    print("  2. Import into OpenFOAM case (constant/triSurface/)")
    print("  3. Run snappyHexMesh with new geometry")
    print()


if __name__ == "__main__":
    main()
