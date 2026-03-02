# UN Headquarters Building Extraction

Demonstration of the CityGML building extraction pipeline applied to the United Nations Headquarters complex on the East River, Manhattan.

## Purpose

This side project demonstrates that the same STL generation pipeline used for the Citicorp Center CFD simulation can extract 3D building geometry for any location in NYC. The UN Headquarters was chosen as a test case because of its distinctive architecture and waterfront setting.

## Files

| File | Description |
|------|-------------|
| `un_headquarters.gml` | Extracted CityGML LOD2 geometry (DA12 district) |
| `un_headquarters.stl` | Triangulated STL mesh for visualization |
| `un_hq_surroundings.stl` | Surrounding buildings within extraction radius |
| `query_un_addresses.py` | Address lookup utility for UN HQ area |

## How to Reproduce

```bash
# Extract UN HQ and surrounding buildings
python ../tools/extract_buildings.py \
    --lat 40.7489 --lon -73.9680 \
    --radius 300 --name un_headquarters

# Or use the hybrid generator for a full CFD-ready output
python ../generators/generate_stl_hybrid.py \
    --center-lat 40.7489 --center-lon -73.9680 \
    --radius 300 --name un_headquarters
```

## Visualizations

- Screenshots: `../../screenshots/un_hq/` (ParaView renderings from multiple angles)
- Presentation: `../../presentations/un_hq_presentation.pptx`
- ParaView script: `../paraview/paraview_un_hq.py`

## Data Source

NYC 3D Building Model (v20v5) — photogrammetric survey from August 2014 aerial imagery.
Source: [georocket/new-york-city-model-enhanced](https://github.com/georocket/new-york-city-model-enhanced)
