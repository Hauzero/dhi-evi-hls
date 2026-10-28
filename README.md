# DHI / EVI from HLS (Google Earth Engine)

Exports Dynamic Habitat Index (DHI) layers, computed from the Enhanced
Vegetation Index (EVI) of NASA HLS (Sentinel-2 + Landsat, 30 m), for polygon
tiles defined in a GeoJSON file.

## Method (per tile and year)
1. HLS S30 + L30 (`NASA/HLS/HLSS30/v002`, `HLSL30/v002`), EVI per image.
2. Fmask: clouds, cloud shadow, adjacent cloud and water masked; snow set to 0.
3. 12 monthly median composites using data from year-1 to year+1.
4. Persistent water (JRC GSW occurrence > 50 %) masked.
5. Gap-fill from neighbouring months, then Savitzky-Golay smoothing (window 5, order 2), clamped to 0-1.
6. DHI components over the 12 months. Output bands (float32, COG GeoTIFF):
   `DHI_sum` (cumulative productivity), `DHI_min` (minimum), `DHI_stdDev` (seasonality), `EVI_median`.

Default years: 2018, 2021, 2024.

## Usage
1. `pip install -r requirements.txt`
2. You need a Google Earth Engine account and your own Cloud project
   ([register](https://earthengine.google.com/)).
3. Run from this folder:

```python
import ee
import dhi_evi

ee.Authenticate()                      # first time only
ee.Initialize(project="your-gee-cloud-project")

tiles = dhi_evi.load_tiles("tiles/proposed_additional_tif_tiles.geojson")
# tiles = dhi_evi.load_tiles(..., tile_ids=["ADD03", "ADD04"])   # subset
dhi_evi.run_exports(tiles, years=[2018, 2021, 2024], drive_folder="DHI_tiles")
```

Exports land in Google Drive under `DHI_tiles/<year>/<tile_id>_DHI_<year>.tif`.
Progress: https://code.earthengine.google.com/tasks

To find and resubmit failed tasks afterwards:

```python
failed = dhi_evi.failed_tasks({"ADD01", "ADD02"}, [2018, 2021, 2024])
by_id = {ft["properties"]["tile_id"]: ft for ft in tiles}
for (tid, yr) in failed:
    dhi_evi.export_tile_year(tid, yr, dhi_evi.build_roi(by_id[tid]), "DHI_tiles")
```

`tiles/proposed_additional_tif_tiles.geojson` holds the 10 additional tiles (ADD01-ADD10).

## Notes
- The largest tiles (ADD01, ADD02, ADD07; up to ~180 x 170 km) may exceed GEE
  export limits. If they fail repeatedly, split the tile into smaller polygons.

## Background
The project concept note is in `docs/`.

## Citation
Code and data: Dr. Leon Hauser, Department of Geography, University of Zurich (UZH).
Please credit the author when using the code or the derived DHI/EVI data.

## Licence
MIT, see `LICENSE`.
