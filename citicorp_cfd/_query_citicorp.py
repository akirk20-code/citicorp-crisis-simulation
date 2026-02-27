"""Query NYC Open Data for the Citicorp Center building footprint."""
import requests, json, math, sys

lat, lon = 40.7585, -73.9702
delta = 0.002
url = "https://data.cityofnewyork.us/resource/5zhs-2jue.geojson"
params = {
    "$where": "within_box(the_geom, %f, %f, %f, %f)" % (lat-delta, lon-delta, lat+delta, lon+delta),
    "$limit": "200"
}
r = requests.get(url, params=params)
data = r.json()
features = data.get("features", [])

citicorp = None
for f in features:
    p = f.get("properties", {})
    if p.get("bin") == "1036474":
        citicorp = f
        break

if not citicorp:
    print("Citicorp not found!")
    sys.exit(1)

p = citicorp["properties"]
geom = citicorp["geometry"]
print("=== Citigroup Center (BIN 1036474) ===")
for k, v in sorted(p.items()):
    print("  %s: %s" % (k, v))

coords = geom["coordinates"]
ring = coords[0][0] if geom["type"] == "MultiPolygon" else coords[0]
print("\nGeometry type: %s" % geom["type"])
print("Vertices: %d" % len(ring))

# Convert to meters relative to centroid
cx = sum(pt[0] for pt in ring) / len(ring)
cy = sum(pt[1] for pt in ring) / len(ring)
m_per_deg_lon = 111320 * math.cos(math.radians(cy))
m_per_deg_lat = 111320

print("\nCentroid: lon=%.7f, lat=%.7f" % (cx, cy))
print("\nFootprint vertices (meters from centroid):")
for i, pt in enumerate(ring):
    x = (pt[0] - cx) * m_per_deg_lon
    y = (pt[1] - cy) * m_per_deg_lat
    print("  [%2d] x=%8.2f m, y=%8.2f m" % (i, x, y))

xs = [(pt[0] - cx) * m_per_deg_lon for pt in ring]
ys = [(pt[1] - cy) * m_per_deg_lat for pt in ring]
print("\nBounding box: X=[%.1f, %.1f] m  Y=[%.1f, %.1f] m" % (min(xs), max(xs), min(ys), max(ys)))
print("Width X: %.1f m  Width Y: %.1f m" % (max(xs) - min(xs), max(ys) - min(ys)))

# Save raw GeoJSON for reference
with open("citicorp_footprint.json", "w") as f:
    json.dump(citicorp, f, indent=2)
print("\nSaved raw GeoJSON to citicorp_footprint.json")
