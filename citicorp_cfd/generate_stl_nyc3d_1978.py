#!/usr/bin/env python3
"""
Generate STL geometry for 1978 historical Citicorp CFD simulation.

Uses NYC Open Data Building Footprints API with construction_year <= 1978
filter. This is the most reliable source for historical data.

Key differences from 2024 version:
  - Only includes buildings confirmed to exist by 1978
  - Excludes 432 Park Ave (2015), Bloomberg Tower (2005), etc.
  - Mean neighbor height: 29.4m (vs 36.0m for 2024)
  - 1,472 buildings (vs 1,665 for 2024)

Citicorp geometry: hand-crafted with 45-deg slanted roof + 45-deg rotation.
All dimensions in meters. Citicorp at origin.
"""

import os, sys, math, struct, json

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

# ─── CONSTANTS ────────────────────────────────────────────────────────────────
CITICORP_LAT = 40.7579
CITICORP_LON = -73.9690

TOWER_W      = 47.85
TOWER_H_TOP  = 278.90
ROOF_START   = 248.00
STILT_H      = 34.75
STILT_W      = 7.32
HALF_T       = TOWER_W / 2
HALF_S       = STILT_W / 2
STILTS       = [(0, -HALF_T), (HALF_T, 0), (0, HALF_T), (-HALF_T, 0)]

M_PER_DEG_LAT = 111320.0
M_PER_DEG_LON = M_PER_DEG_LAT * math.cos(math.radians(CITICORP_LAT))
FT_TO_M       = 0.3048
MIN_HEIGHT_M  = 5.0
HISTORICAL_YEAR = 1978

BB_LAT_MIN, BB_LAT_MAX = 40.7507, 40.7651
BB_LON_MIN, BB_LON_MAX = -73.9762, -73.9618

NYC_API = "https://data.cityofnewyork.us/resource/5zhs-2jue.geojson"
CITICORP_BINS = {"1035879", "1087931"}

# ─── GEOMETRY HELPERS ─────────────────────────────────────────────────────────
def wgs84_to_local(lon, lat):
    return (lon - CITICORP_LON) * M_PER_DEG_LON, (lat - CITICORP_LAT) * M_PER_DEG_LAT

def rotate45(x, y):
    c = s = math.sqrt(0.5)
    return c*x - s*y, s*x + c*y

def poly_area_2d(v):
    n = len(v)
    a = 0.0
    for i in range(n):
        j = (i+1)%n
        a += v[i][0]*v[j][1] - v[j][0]*v[i][1]
    return a/2.0

def ear_clip(verts):
    pts = list(verts)
    if poly_area_2d(pts) < 0:
        pts.reverse()
    tris = []
    while len(pts) > 3:
        n = len(pts)
        found = False
        for i in range(n):
            a,b,c = pts[(i-1)%n], pts[i], pts[(i+1)%n]
            cross = (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
            if cross <= 0:
                continue
            ok = True
            for j in range(n):
                if j in ((i-1)%n, i, (i+1)%n):
                    continue
                p = pts[j]
                d1=(b[0]-a[0])*(p[1]-a[1])-(b[1]-a[1])*(p[0]-a[0])
                d2=(c[0]-b[0])*(p[1]-b[1])-(c[1]-b[1])*(p[0]-b[0])
                d3=(a[0]-c[0])*(p[1]-c[1])-(a[1]-c[1])*(p[0]-c[0])
                if (d1>=0 and d2>=0 and d3>=0) or (d1<=0 and d2<=0 and d3<=0):
                    ok = False; break
            if ok:
                tris.append((a,b,c)); pts.pop(i); found = True; break
        if not found:
            break
    if len(pts)==3:
        tris.append(tuple(pts))
    return tris

def tri_normal(p0,p1,p2):
    ax,ay,az = p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2]
    bx,by,bz = p2[0]-p0[0],p2[1]-p0[1],p2[2]-p0[2]
    nx,ny,nz = ay*bz-az*by, az*bx-ax*bz, ax*by-ay*bx
    m = math.sqrt(nx*nx+ny*ny+nz*nz)
    return (nx/m,ny/m,nz/m) if m>1e-12 else (0.,0.,1.)

def write_stl(path, name, tris):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path,'w') as f:
        f.write(f"solid {name}\n")
        for t in tris:
            n = tri_normal(*t)
            f.write(f"  facet normal {n[0]:.6f} {n[1]:.6f} {n[2]:.6f}\n    outer loop\n")
            for p in t:
                f.write(f"      vertex {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
            f.write("    endloop\n  endfacet\n")
        f.write(f"endsolid {name}\n")
    print(f"  Wrote {len(tris):6d} triangles -> {os.path.basename(path)}")
    return len(tris)

def extrude(verts2d, zbot, ztop):
    tris = []
    v = list(verts2d)
    if poly_area_2d(v) < 0: v.reverse()
    n = len(v)
    for i in range(n):
        a,b = v[i],v[(i+1)%n]
        p00=(a[0],a[1],zbot); p10=(b[0],b[1],zbot)
        p11=(b[0],b[1],ztop); p01=(a[0],a[1],ztop)
        tris += [(p00,p10,p11),(p00,p11,p01)]
    for t in ear_clip(v):
        tris.append(((t[0][0],t[0][1],ztop),(t[1][0],t[1][1],ztop),(t[2][0],t[2][1],ztop)))
    for t in ear_clip(v):
        tris.append(((t[0][0],t[0][1],zbot),(t[2][0],t[2][1],zbot),(t[1][0],t[1][1],zbot)))
    return tris

# ─── CITICORP GEOMETRY ────────────────────────────────────────────────────────
def make_tower():
    hw = HALF_T
    r2 = math.sqrt(2.0)
    corners = [rotate45(x,y) for x,y in [(-hw,-hw),(hw,-hw),(hw,hw),(-hw,hw)]]
    y_s = -hw*r2; y_n = hw*r2; yr = y_n - y_s

    def rz(cy): return ROOF_START + (cy - y_s)/yr * (TOWER_H_TOP - ROOF_START)

    tris = []
    zb, zt = STILT_H, ROOF_START

    # Walls (box section)
    for i in range(4):
        c0,c1 = corners[i],corners[(i+1)%4]
        p00=(c0[0],c0[1],zb); p10=(c1[0],c1[1],zb)
        p11=(c1[0],c1[1],zt); p01=(c0[0],c0[1],zt)
        tris += [(p00,p10,p11),(p00,p11,p01)]

    # Bottom cap
    base2d = [(c[0],c[1]) for c in corners]
    for t in ear_clip(base2d):
        tris.append(((t[0][0],t[0][1],zb),(t[2][0],t[2][1],zb),(t[1][0],t[1][1],zb)))

    # Slanted roof (planar)
    rp = [(c[0],c[1],rz(c[1])) for c in corners]
    tris += [(rp[0],rp[1],rp[2]),(rp[0],rp[2],rp[3])]

    # Transition triangles between box top and slanted roof
    for i in range(4):
        c0,c1 = corners[i],corners[(i+1)%4]
        rb=(c0[0],c0[1],zt); re=(c1[0],c1[1],zt)
        rt0=(c0[0],c0[1],rz(c0[1])); rt1=(c1[0],c1[1],rz(c1[1]))
        if abs(rt0[2]-zt) > 0.01 or abs(rt1[2]-zt) > 0.01:
            tris.append((rb,re,rt1))
            if abs(rt0[2]-zt) > 0.01:
                tris.append((rb,rt1,rt0))
    return tris

def make_stilts():
    tris = []
    for cx,cy in STILTS:
        rx,ry = rotate45(cx,cy)
        hw = HALF_S
        c = [(rx-hw,ry-hw),(rx+hw,ry-hw),(rx+hw,ry+hw),(rx-hw,ry+hw)]
        zb,zt = 0.0, STILT_H
        for i in range(4):
            a,b = c[i],c[(i+1)%4]
            tris += [((a[0],a[1],zb),(b[0],b[1],zb),(b[0],b[1],zt)),
                     ((a[0],a[1],zb),(b[0],b[1],zt),(a[0],a[1],zt))]
        tris += [((c[0][0],c[0][1],zt),(c[2][0],c[2][1],zt),(c[1][0],c[1][1],zt)),
                 ((c[0][0],c[0][1],zt),(c[3][0],c[3][1],zt),(c[2][0],c[2][1],zt))]
    return tris

# ─── NYC API WITH 1978 FILTER ─────────────────────────────────────────────────
def fetch_1978():
    print(f"  Querying NYC Building Footprints API (construction_year <= {HISTORICAL_YEAR})...")
    where = (f"within_box(the_geom, {BB_LAT_MAX}, {BB_LON_MIN}, {BB_LAT_MIN}, {BB_LON_MAX})"
             f" AND (construction_year IS NULL OR construction_year <= {HISTORICAL_YEAR})")
    r = requests.get(NYC_API, params={"$where": where, "$limit": 10000}, timeout=60)
    r.raise_for_status()
    features = r.json().get("features", [])
    print(f"  Received {len(features)} buildings (pre-{HISTORICAL_YEAR})")
    return features

def parse_feature(feat):
    props = feat.get("properties", {})
    geom  = feat.get("geometry", {})
    hft   = props.get("height_roof")
    if not hft: return None
    hm = float(hft)*FT_TO_M
    if hm < MIN_HEIGHT_M: return None
    gtype = geom.get("type","")
    coords= geom.get("coordinates",[])
    ring  = None
    if gtype=="Polygon" and coords: ring=coords[0]
    elif gtype=="MultiPolygon" and coords:
        best,ba=None,0
        for p in coords:
            if p and p[0]:
                a=abs(poly_area_2d([(c[0],c[1]) for c in p[0]]))
                if a>ba: ba,best=a,p[0]
        ring=best
    if not ring or len(ring)<4: return None
    verts=[wgs84_to_local(c[0],c[1]) for c in ring]
    if math.hypot(verts[0][0]-verts[-1][0],verts[0][1]-verts[-1][1])<0.01:
        verts.pop()
    if len(verts)<3 or abs(poly_area_2d(verts))<10: return None
    return verts, hm

def is_citicorp(feat):
    props = feat.get("properties",{})
    if str(props.get("bin","")) in CITICORP_BINS: return True
    geom  = feat.get("geometry",{})
    gtype = geom.get("type","")
    coords= geom.get("coordinates",[])
    ring  = None
    if gtype=="Polygon" and coords: ring=coords[0]
    elif gtype=="MultiPolygon" and coords and coords[0]: ring=coords[0][0]
    if ring:
        v = [wgs84_to_local(c[0],c[1]) for c in ring]
        if v:
            cx=sum(p[0] for p in v)/len(v)
            cy=sum(p[1] for p in v)/len(v)
            if math.hypot(cx,cy)<35: return True
    return False

# ─── MAIN ─────────────────────────────────────────────────────────────────────
def main():
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "constant","triSurface")
    os.makedirs(out_dir, exist_ok=True)

    print("\nCiticorp 1978 Historical — NYC Open Data Source")
    print("="*60)

    print("\nBuilding Citicorp tower (slanted roof, 45-deg rotation)...")
    nt = write_stl(os.path.join(out_dir,"citicorp_tower.stl"), "citicorp_tower", make_tower())

    print("\nBuilding Citicorp stilts...")
    ns = write_stl(os.path.join(out_dir,"citicorp_stilts.stl"), "citicorp_stilts", make_stilts())

    print(f"\nFetching pre-{HISTORICAL_YEAR} surroundings...")
    surr_tris=[]; count=0; heights=[]

    if HAS_REQUESTS:
        try:
            features = fetch_1978()
            for feat in features:
                if is_citicorp(feat): continue
                p = parse_feature(feat)
                if not p: continue
                verts, hm = p
                surr_tris.extend(extrude(verts, 0.0, hm))
                count += 1
                heights.append(hm)
        except Exception as e:
            print(f"  API error: {e}")
    else:
        print("  requests not installed (pip install requests)")

    nb = write_stl(os.path.join(out_dir,"surroundings.stl"),
                   f"surroundings_nyc3d_{HISTORICAL_YEAR}", surr_tris)

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"1978 HISTORICAL BUILD SUMMARY (NYC Open Data source)")
    print(f"{'='*60}")
    print(f"Tower triangles:         {nt}")
    print(f"Stilt triangles:         {ns}")
    print(f"Surrounding triangles:   {nb}")
    print(f"Buildings included:      {count}")
    if heights:
        print(f"Min/Max height:          {min(heights):.1f}m / {max(heights):.1f}m")
        print(f"Mean height:             {sum(heights)/len(heights):.1f}m")
        print(f"Tallest known (1978):    Chrysler ~288.6m")
        print(f"Citicorp rank (1978):    2nd (after Chrysler)")

    print(f"\nKey 1978 vs 2024 differences:")
    print(f"  Buildings: {count} (vs ~1472 total pre-1978)")
    print(f"  Excludes:  432 Park Ave (2015), Bloomberg Tower (2005), etc.")
    print(f"\nOutput: {out_dir}")
    print(f"\nNext: copy STLs to ~/citicorp_cfd_nyc3d_1978/constant/triSurface/")

if __name__ == "__main__":
    main()
