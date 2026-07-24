#!/usr/bin/env python3
"""
Download a satellite texture for an OpenPDAC/OpenFOAM DEM domain.

Key rule:
    The source CRS is mandatory. The script never guesses it from the DEM and
    never uses a silent fallback. Pass the correct EPSG code explicitly with
    --epsg.

Important note:
    Web tile providers used by contextily serve imagery in Web Mercator
    (EPSG:3857). Therefore, the script must internally transform the DEM bounds
    to EPSG:3857 only to request the imagery tiles. The final texture is then
    resampled back onto the DEM grid defined by the user-provided EPSG code.

Outputs:
    - A JPEG texture image matching the DEM extent and aspect ratio.
    - The output texture resolution can be increased with --scale, or set
      explicitly with --width and/or --height.
    - Optionally, a GeoTIFF debug file in the user-provided CRS, useful for GIS
      verification.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Tuple

import contextily as cx
import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.transform import from_bounds
from rasterio.warp import reproject, transform_bounds
from PIL import Image


WEB_MERCATOR_CRS = CRS.from_epsg(3857)


PROVIDERS = {
    "esri": cx.providers.Esri.WorldImagery,
    "osm": cx.providers.OpenStreetMap.Mapnik,
    "opentopo": cx.providers.OpenTopoMap,
}


def parse_epsg(value: str) -> CRS:
    """Parse an explicit EPSG code passed by the user."""
    text = value.strip().upper()
    if text.startswith("EPSG:"):
        code = text.split(":", 1)[1]
    else:
        code = text

    if not code.isdigit():
        raise argparse.ArgumentTypeError(
            "The EPSG code must be written as an integer, for example 32633, "
            "or as EPSG:32633."
        )

    return CRS.from_epsg(int(code))


def parse_zoom(value: str):
    """Accept either 'auto' or an integer zoom level."""
    if value.lower() == "auto":
        return "auto"
    try:
        zoom = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Zoom must be 'auto' or an integer.") from exc
    if zoom < 0 or zoom > 23:
        raise argparse.ArgumentTypeError("Zoom must be between 0 and 23.")
    return zoom


def get_provider(name: str):
    """Return a contextily tile provider from a short provider name."""
    key = name.lower()
    if key not in PROVIDERS:
        valid = ", ".join(sorted(PROVIDERS))
        raise argparse.ArgumentTypeError(f"Unknown provider '{name}'. Valid providers: {valid}.")
    return PROVIDERS[key]


def read_dem_metadata(dem_path: Path, source_crs: CRS) -> Tuple[rasterio.coords.BoundingBox, object, int, int]:
    """Read DEM bounds, affine transform, width and height."""
    if not dem_path.exists():
        raise FileNotFoundError(f"DEM file not found: {dem_path}")

    print(f"1. Reading DEM: {dem_path}")
    with rasterio.open(dem_path) as src:
        bounds = src.bounds
        transform = src.transform
        width = src.width
        height = src.height
        embedded_crs = src.crs

    print(f"   - User-provided source CRS: {source_crs.to_string()}")
    if embedded_crs:
        print(f"   - CRS embedded in DEM, ignored by design: {embedded_crs.to_string()}")
    else:
        print("   - DEM has no embedded CRS. This is OK because --epsg is mandatory.")
    print(f"   - DEM size: {width} x {height} pixels")
    print(
        "   - DEM bounds in user-provided CRS: "
        f"({bounds.left:.6f}, {bounds.bottom:.6f}) -> ({bounds.right:.6f}, {bounds.top:.6f})"
    )

    return bounds, transform, width, height


def download_tile_mosaic(bounds, source_crs: CRS, provider, zoom):
    """Download a Web Mercator tile mosaic covering the DEM bounds."""
    print("\n2. Preparing tile request")
    print("   - The DEM CRS remains the user-provided CRS.")
    print("   - Only the tile request is transformed internally to EPSG:3857, as required by web tile servers.")

    west_3857, south_3857, east_3857, north_3857 = transform_bounds(
        source_crs,
        WEB_MERCATOR_CRS,
        bounds.left,
        bounds.bottom,
        bounds.right,
        bounds.top,
        densify_pts=21,
    )

    print(
        "   - Tile request bounds in EPSG:3857: "
        f"({west_3857:.3f}, {south_3857:.3f}) -> ({east_3857:.3f}, {north_3857:.3f})"
    )
    print(f"   - Provider: {provider.name}")
    print(f"   - Zoom: {zoom}")

    img, extent_3857 = cx.bounds2img(
        west_3857,
        south_3857,
        east_3857,
        north_3857,
        ll=False,
        source=provider,
        zoom=zoom,
    )

    print("   - Tile download completed")
    print(f"   - Downloaded mosaic size: {img.shape[1]} x {img.shape[0]} pixels")
    print(
        "   - Downloaded mosaic extent in EPSG:3857: "
        f"left={extent_3857[0]:.3f}, right={extent_3857[1]:.3f}, "
        f"bottom={extent_3857[2]:.3f}, top={extent_3857[3]:.3f}"
    )

    return img, extent_3857


def resample_mosaic_to_target_grid(
    img: np.ndarray,
    extent_3857,
    target_transform,
    target_width: int,
    target_height: int,
    source_crs: CRS,
) -> np.ndarray:
    """Resample the downloaded EPSG:3857 mosaic onto the target texture grid."""
    print("\n4. Resampling tile mosaic onto the target texture grid")
    print("   - Target grid CRS: user-provided DEM CRS")
    print("   - Target grid extent and aspect ratio: DEM extent")

    left_3857, right_3857, bottom_3857, top_3857 = extent_3857
    src_height, src_width = img.shape[:2]
    src_transform = from_bounds(
        left_3857,
        bottom_3857,
        right_3857,
        top_3857,
        src_width,
        src_height,
    )

    if img.ndim != 3:
        raise ValueError("The downloaded image is expected to have RGB or RGBA bands.")

    # Keep RGB only. JPEG does not support alpha in the usual way.
    src_rgb = img[:, :, :3]
    dst_rgb = np.zeros((target_height, target_width, 3), dtype=np.uint8)

    for band in range(3):
        reproject(
            source=src_rgb[:, :, band],
            destination=dst_rgb[:, :, band],
            src_transform=src_transform,
            src_crs=WEB_MERCATOR_CRS,
            dst_transform=target_transform,
            dst_crs=source_crs,
            resampling=Resampling.bilinear,
            dst_nodata=0,
        )

    print(f"   - Final texture size: {target_width} x {target_height} pixels")
    return dst_rgb


def save_jpeg(texture_rgb: np.ndarray, output_path: Path, quality: int):
    """Save the final texture as JPEG."""
    print(f"\n5. Saving JPEG texture: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(texture_rgb).save(output_path, format="JPEG", quality=quality, subsampling=0)
    print("   - JPEG texture saved")


def save_debug_geotiff(texture_rgb: np.ndarray, output_path: Path, transform, crs: CRS):
    """Save a georeferenced debug GeoTIFF in the user-provided CRS."""
    print(f"\n6. Saving debug GeoTIFF: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    height, width = texture_rgb.shape[:2]
    profile = {
        "driver": "GTiff",
        "height": height,
        "width": width,
        "count": 3,
        "dtype": "uint8",
        "crs": crs,
        "transform": transform,
        "compress": "deflate",
        "photometric": "RGB",
    }

    with rasterio.open(output_path, "w", **profile) as dst:
        for band in range(3):
            dst.write(texture_rgb[:, :, band], band + 1)

    print("   - Debug GeoTIFF saved")



def compute_target_grid(bounds, dem_width: int, dem_height: int, scale: float, width, height):
    """Compute the output texture grid while preserving the DEM extent."""
    if scale <= 0:
        raise ValueError("--scale must be greater than zero.")

    aspect = dem_width / dem_height

    if width is not None and width <= 0:
        raise ValueError("--width must be greater than zero.")
    if height is not None and height <= 0:
        raise ValueError("--height must be greater than zero.")

    if width is not None and height is not None:
        target_width = int(width)
        target_height = int(height)
        print("   - Output resolution mode: explicit --width and --height")
    elif width is not None:
        target_width = int(width)
        target_height = max(1, int(round(target_width / aspect)))
        print("   - Output resolution mode: explicit --width, height from DEM aspect ratio")
    elif height is not None:
        target_height = int(height)
        target_width = max(1, int(round(target_height * aspect)))
        print("   - Output resolution mode: explicit --height, width from DEM aspect ratio")
    else:
        target_width = max(1, int(round(dem_width * scale)))
        target_height = max(1, int(round(dem_height * scale)))
        print(f"   - Output resolution mode: DEM grid multiplied by --scale={scale:g}")

    target_transform = from_bounds(
        bounds.left,
        bounds.bottom,
        bounds.right,
        bounds.top,
        target_width,
        target_height,
    )

    return target_transform, target_width, target_height

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download a satellite texture exactly aligned to a DEM extent. The DEM EPSG code is mandatory.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--epsg",
        required=True,
        type=parse_epsg,
        help="Mandatory CRS of the DEM, for example 32633 or EPSG:32633.",
    )
    parser.add_argument(
        "--dem",
        default="./constant/DEM/DEMcropped.asc",
        help="Input DEM file path.",
    )
    parser.add_argument(
        "--output",
        default="./constant/DEM/texture.jpg",
        help="Output JPEG texture path.",
    )
    parser.add_argument(
        "--debug-geotiff",
        default="./constant/DEM/texture_debug_user_crs.tif",
        help="Output debug GeoTIFF path. Use --no-debug-geotiff to disable it.",
    )
    parser.add_argument(
        "--no-debug-geotiff",
        action="store_true",
        help="Do not write the georeferenced debug GeoTIFF.",
    )
    parser.add_argument(
        "--provider",
        default="esri",
        type=get_provider,
        help="Tile provider short name: esri, osm, or opentopo.",
    )
    parser.add_argument(
        "--zoom",
        default="auto",
        type=parse_zoom,
        help="Tile zoom level: 'auto' or an integer.",
    )
    parser.add_argument(
        "--scale",
        default=1.0,
        type=float,
        help=(
            "Multiplier for the DEM pixel grid used as final texture resolution. "
            "For example, --scale 4 writes a texture with 4 times the DEM width "
            "and 4 times the DEM height. Ignored if --width or --height is used."
        ),
    )
    parser.add_argument(
        "--width",
        default=None,
        type=int,
        help=(
            "Explicit output texture width in pixels. If --height is omitted, "
            "height is computed from the DEM aspect ratio."
        ),
    )
    parser.add_argument(
        "--height",
        default=None,
        type=int,
        help=(
            "Explicit output texture height in pixels. If --width is omitted, "
            "width is computed from the DEM aspect ratio."
        ),
    )
    parser.add_argument(
        "--jpeg-quality",
        default=95,
        type=int,
        help="JPEG quality, from 1 to 100.",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.jpeg_quality < 1 or args.jpeg_quality > 100:
        parser.error("--jpeg-quality must be between 1 and 100.")

    dem_path = Path(args.dem)
    output_path = Path(args.output)
    debug_geotiff_path = Path(args.debug_geotiff)
    source_crs = args.epsg

    print("Satellite texture downloader for DEM-based rendering")
    print("====================================================")
    print("CRS policy: the DEM CRS must be passed explicitly with --epsg.")
    print("No CRS fallback is used. No CRS guessed from the DEM is used.")

    bounds, dem_transform, dem_width, dem_height = read_dem_metadata(dem_path, source_crs)
    print("\n3. Choosing output texture resolution")
    target_transform, target_width, target_height = compute_target_grid(
        bounds,
        dem_width,
        dem_height,
        args.scale,
        args.width,
        args.height,
    )
    print(f"   - DEM size: {dem_width} x {dem_height} pixels")
    print(f"   - Output texture size: {target_width} x {target_height} pixels")

    img, extent_3857 = download_tile_mosaic(bounds, source_crs, args.provider, args.zoom)
    texture_rgb = resample_mosaic_to_target_grid(
        img,
        extent_3857,
        target_transform,
        target_width,
        target_height,
        source_crs,
    )
    save_jpeg(texture_rgb, output_path, args.jpeg_quality)

    if not args.no_debug_geotiff:
        save_debug_geotiff(texture_rgb, debug_geotiff_path, target_transform, source_crs)

    print("\nDone.")
    print(f"Texture image: {output_path}")
    if not args.no_debug_geotiff:
        print(f"Debug GeoTIFF: {debug_geotiff_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
