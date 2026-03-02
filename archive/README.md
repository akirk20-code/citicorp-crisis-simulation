# Archive — Development Iterations

This directory preserves superseded scripts and files that show the project's development evolution. These files are no longer actively used but are retained for reference.

## Generators

Superseded STL generation scripts. The project explored multiple data sources and approaches before settling on the NYC Open Data API (`generators/generate_stl_nyc3d.py`).

| Script | Approach | Why Superseded |
|--------|----------|----------------|
| `generate_stl_osm.py` | OpenStreetMap Overpass API | Variable height accuracy, incorrect stilt positions |
| `generate_stl_osm_1978.py` | OSM with 1978 filter | Same issues as OSM variant |
| `generate_stl_lod2.py` | TUM CityGML LOD2 | Complex setup, falls back to API anyway |
| `generate_stl_lod2_1978.py` | CityGML with 1978 filter | Superseded by `generate_stl_hybrid.py --year 1978` |
| `generate_stl_nyc3d_1978.py` | NYC API with 1978 filter | Superseded by `generate_stl_hybrid.py --year 1978` |

## Screenshots

Earlier versions of comparison screenshots, superseded by v3 iterations in `screenshots/cityscape/`.

## Presentations

Earlier versions of PowerPoint files, superseded by v2 in `presentations/`.

---

**Note:** The full git history is preserved. To see any file at any point in time:
```bash
git log --follow -- archive/generators/generate_stl_osm.py
git show v0.1-pre-reorganization:citicorp_cfd/generate_stl_osm.py
```
