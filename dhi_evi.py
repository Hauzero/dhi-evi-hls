"""
DHI / EVI export from HLS (Sentinel-2 + Landsat) via Google Earth Engine.

For each tile polygon and target year, builds 12 monthly median EVI composites
(+/- 1 year of data), gap-fills and Savitzky-Golay smooths them, and exports
four bands to Google Drive as a Cloud-Optimized GeoTIFF:
    DHI_sum, DHI_min, DHI_stdDev, EVI_median
"""
import json
import time

import ee
from scipy import signal

# =============================================================================
# DEFAULT CONFIGURATION (override via keyword arguments)
# =============================================================================
TARGET_YEARS = [2018, 2021, 2024]
EXPORT_SCALE = 30          # metres
EXPORT_FLOAT32 = True      # halves file size; values unchanged to ~7 digits
MAX_VERTICES = 3000        # simplify rings above this before sending to GEE
SIMPLIFY_TOLERANCE = 0.01  # degrees (~1 km)
BATCH_SIZE = 20            # tiles per submission batch (rate limiting)


# =============================================================================
# TILE LOADING / GEOMETRY
# =============================================================================
def load_tiles(path, tile_ids=None):
    """Load features from a GeoJSON file, optionally keeping only some tile_ids."""
    with open(path, encoding="utf-8") as f:
        features = json.load(f)["features"]
    if tile_ids:
        wanted = set(tile_ids)
        features = [ft for ft in features if ft["properties"].get("tile_id") in wanted]
    return features


def _simplify_ring(coords, tolerance=SIMPLIFY_TOLERANCE):
    try:
        from shapely.geometry import LinearRing
        simplified = LinearRing(coords).simplify(tolerance, preserve_topology=False)
        result = list(simplified.coords)
    except ImportError:
        step = max(1, len(coords) // MAX_VERTICES)
        result = coords[::step]
        if result[0] != result[-1]:
            result.append(result[0])
    return [[c[0], c[1]] for c in result]


def build_roi(feature):
    """Return an ee.Geometry.Polygon for a GeoJSON Polygon feature."""
    ring = feature["geometry"]["coordinates"][0]
    if len(ring) > MAX_VERTICES:
        ring = _simplify_ring(ring)
    else:
        ring = [[c[0], c[1]] for c in ring]
    return ee.Geometry.Polygon([ring])


# =============================================================================
# DATA PROCESSING
# =============================================================================
def add_evi_hls(image):
    """Enhanced Vegetation Index for HLS data."""
    evi = image.expression(
        "2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))",
        {"NIR": image.select("NIR"), "RED": image.select("Red"), "BLUE": image.select("Blue")},
    ).rename("EVI")
    return image.addBands(evi)


def mask_hls_pixels(image):
    """Mask cloud/shadow/adjacent/water via Fmask; set snow to 0. Returns EVI."""
    fmask = image.select("Fmask")
    cloud_shadow, clouds, adjacent = 1 << 3, 1 << 1, 1 << 2
    snow, water = 1 << 4, 1 << 5

    cloud_mask = (fmask.bitwiseAnd(cloud_shadow).eq(0)
                  .And(fmask.bitwiseAnd(clouds).eq(0))
                  .And(fmask.bitwiseAnd(adjacent).eq(0)))
    snow_mask = fmask.bitwiseAnd(snow).eq(snow)
    water_mask = fmask.bitwiseAnd(water).eq(water)

    evi = image.select("EVI").where(snow_mask, 0)
    evi = evi.updateMask(cloud_mask).updateMask(water_mask.Not())
    return evi.rename("EVI")


def create_jrc_water_mask():
    """Persistent water from JRC Global Surface Water (>50% occurrence)."""
    occurrence = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("occurrence")
    return occurrence.gt(50).unmask(ee.Image(0)).rename("persistent_water")


def create_monthly_composite(month_num, hls_combined, persistent_water, target_year):
    """Median monthly composite using +/-1 year of HLS data."""
    images = ee.ImageCollection([])
    for offset in (-1, 0, 1):
        start = ee.Date.fromYMD(target_year + offset, month_num, 1)
        images = images.merge(
            hls_combined.filterDate(start, start.advance(1, "month")).select("EVI"))
    composite = images.median().updateMask(persistent_water.Not())
    return (composite
            .set("month", month_num)
            .set("year", target_year)
            .set("hls_images", images.size())
            .set("system:time_start", ee.Date.fromYMD(target_year, month_num, 1).millis()))


def _dhi_reducer():
    """sum / min / stdDev / median in a single pass."""
    return (ee.Reducer.sum()
            .combine(ee.Reducer.min(), sharedInputs=True)
            .combine(ee.Reducer.stdDev(), sharedInputs=True)
            .combine(ee.Reducer.median(), sharedInputs=True))


def compute_dhi(target_year, roi, float32=EXPORT_FLOAT32):
    """Return the 4-band DHI ee.Image for one polygon and year."""
    start_date = ee.Date.fromYMD(target_year - 1, 1, 1)
    end_date = ee.Date.fromYMD(target_year + 1, 12, 31)

    sentinel2 = (ee.ImageCollection("NASA/HLS/HLSS30/v002")
                 .filterBounds(roi).filter(ee.Filter.date(start_date, end_date))
                 .select(["B2", "B4", "B8A", "Fmask"], ["Blue", "Red", "NIR", "Fmask"])
                 .map(add_evi_hls).map(mask_hls_pixels))
    landsat = (ee.ImageCollection("NASA/HLS/HLSL30/v002")
               .filterBounds(roi).filter(ee.Filter.date(start_date, end_date))
               .select(["B2", "B4", "B5", "Fmask"], ["Blue", "Red", "NIR", "Fmask"])
               .map(add_evi_hls).map(mask_hls_pixels))
    hls_combined = sentinel2.merge(landsat)
    persistent_water = create_jrc_water_mask()

    monthly = ee.ImageCollection([
        create_monthly_composite(m, hls_combined, persistent_water, target_year)
        for m in range(1, 13)])

    monthly_list = monthly.toList(12)
    all_valid = monthly.qualityMosaic("EVI")
    originals = [ee.Image(monthly_list.get(i)) for i in range(12)]

    # Gap-fill: average of neighbours at increasing radius, then best-EVI mosaic
    gap_filled = []
    for i in range(12):
        filled = originals[i]
        for radius in range(1, 7):
            if radius < 6:
                prev = originals[(i - radius) % 12]
                nxt = originals[(i + radius) % 12]
                filled = filled.unmask(prev.add(nxt).divide(2))
            else:
                filled = filled.unmask(all_valid)
        gap_filled.append(filled.updateMask(persistent_water.Not()).set("month", i + 1))

    # Savitzky-Golay smoothing (window 5, order 2), cyclic over the year
    sg_coeffs = signal.savgol_coeffs(5, 2, deriv=0)
    padded = gap_filled[-2:] + gap_filled + gap_filled[:2]
    smoothed = []
    for i in range(12):
        smooth = ee.Image(0)
        for img, coeff in zip(padded[i:i + 5], sg_coeffs):
            smooth = smooth.add(img.multiply(float(coeff)))
        smoothed.append(smooth.set("month", i + 1))

    cleaned = ee.ImageCollection([
        img.rename("EVI").clamp(0, 1).updateMask(img.mask())
        .updateMask(persistent_water.Not()).set("month", i + 1)
        for i, img in enumerate(smoothed)])

    dhi = (cleaned.select("EVI").reduce(_dhi_reducer())
           .select(["EVI_sum", "EVI_min", "EVI_stdDev", "EVI_median"],
                   ["DHI_sum", "DHI_min", "DHI_stdDev", "EVI_median"])
           .updateMask(persistent_water.Not())
           .clip(roi))
    return dhi.toFloat() if float32 else dhi


def export_tile_year(tile_id, year, roi, drive_folder="DHI_tiles", scale=EXPORT_SCALE):
    """Start a GEE export task to Drive; returns the output name."""
    name = f"{tile_id}_DHI_{year}"
    ee.batch.Export.image.toDrive(
        image=compute_dhi(year, roi),
        description=f"{name}_GeoTIFF",
        folder=f"{drive_folder}/{year}",
        fileNamePrefix=name,
        region=roi,
        scale=scale,
        maxPixels=1e10,
        fileFormat="GeoTIFF",
        formatOptions={"cloudOptimized": True},
    ).start()
    return name


def run_exports(tiles, years=TARGET_YEARS, drive_folder="DHI_tiles", batch_size=BATCH_SIZE):
    """Submit exports for all tiles x years. Returns (submitted, errors)."""
    submitted, errors = [], []
    for start in range(0, len(tiles), batch_size):
        for ft in tiles[start:start + batch_size]:
            tile_id = ft["properties"]["tile_id"]
            roi = build_roi(ft)
            for year in years:
                try:
                    submitted.append(export_tile_year(tile_id, year, roi, drive_folder))
                    print(f"  OK   {tile_id} / {year}")
                except Exception as e:  # noqa: BLE001
                    errors.append((tile_id, year, str(e)))
                    print(f"  FAIL {tile_id} / {year}: {e}")
        if start + batch_size < len(tiles):
            time.sleep(2)
    print(f"\nSubmitted {len(submitted)} task(s), {len(errors)} error(s).")
    print("Monitor: https://code.earthengine.google.com/tasks")
    return submitted, errors


def failed_tasks(tile_ids, years):
    """List (tile_id, year) whose latest export FAILED and never COMPLETED."""
    import re
    pat = re.compile(r"^(.+)_DHI_(\d{4})_GeoTIFF$")
    done, failed = set(), {}
    for t in ee.batch.Task.list():
        s = t.status()
        m = pat.match(s.get("description", ""))
        if not m or m.group(1) not in tile_ids or int(m.group(2)) not in years:
            continue
        key = (m.group(1), int(m.group(2)))
        if s["state"] == "COMPLETED":
            done.add(key)
        elif s["state"] == "FAILED":
            failed[key] = s.get("error_message", "")
    return {k: v for k, v in failed.items() if k not in done}
