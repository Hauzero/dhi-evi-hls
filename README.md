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
   ([register](https://earthengine.google.com/)). Set `EE_PROJECT` in the notebook.
3. Open `dhi_evi_export.ipynb`, authenticate on first run (`ee.Authenticate()`), run the cells.
   Exports land in Google Drive under `DHI_tiles/<year>/<tile_id>_DHI_<year>.tif`.

`tiles/proposed_additional_tif_tiles.geojson` holds the 10 additional tiles
(ADD01-ADD10). Pass `tile_ids=[...]` to `load_tiles` to run a subset.

## Notes
- The largest tiles (ADD01, ADD02, ADD07; up to ~180 x 170 km) may exceed GEE
  export limits. The notebook has a cell that lists failed tasks; if they
  recur, split the tile into smaller polygons.
- The tiles are exported as the polygon region, not the bounding box.

## Citation
Code and data: Dr. Leon Hauser, Department of Geography, University of Zurich (UZH).
Please credit the author when using the code or the derived DHI/EVI data.

## Licence
MIT, see `LICENSE`.
