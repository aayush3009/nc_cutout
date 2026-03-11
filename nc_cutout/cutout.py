#!/usr/bin/env python3
"""JWST NIRCam cutout creator for JADES survey data.

Creates postage-stamp cutouts from NIRCam imaging FITS files for a list of
sources defined in a CSV catalogue. Supports GOODS-S and GOODS-N fields.
"""

import argparse
import glob
import logging
from pathlib import Path

import astropy.units as u
import pandas as pd
import yaml
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Default configuration
# ---------------------------------------------------------------------------

DEFAULT_CONFIG = {
    "catalogue": {
        "id_column": "ID_1",
        "ra_column": "RA",
        "dec_column": "DEC",
        "survey_column": "SURVEY",
        "goods_s_surveys": [
            "goods-s-deephst",
            "goods-s-deepjwst",
            "goods-s-ultradeep",
            "goods-s-mediumhst",
            "goods-s-mediumjwst",
        ],
        "goods_n_surveys": [
            "goods-n-mediumhst",
            "goods-n-mediumjwst",
        ],
    },
    "images": {
        "sci_extension": "SCI",
        # Index in the underscore-split image filename that gives the filter name.
        # For JADES filenames like ``field_FILTER_...fits`` this is 1.
        "filter_field_index": 1,
    },
    "output": {
        "dir": "./cutouts",
    },
    "cutout": {
        "size_arcsec": 2.0,
    },
}


# ---------------------------------------------------------------------------
# Configuration helpers
# ---------------------------------------------------------------------------

def load_config(config_path=None):
    """Load configuration, merging a user YAML file over the built-in defaults.

    Parameters
    ----------
    config_path : str or Path, optional
        Path to a YAML configuration file.  If ``None`` the default
        configuration is returned unchanged.

    Returns
    -------
    dict
        Merged configuration dictionary.
    """
    import copy
    config = copy.deepcopy(DEFAULT_CONFIG)
    if config_path is not None:
        with open(config_path) as f:
            user_config = yaml.safe_load(f)
        for section, values in user_config.items():
            if section in config and isinstance(config[section], dict):
                config[section].update(values)
            else:
                config[section] = values
    return config


# ---------------------------------------------------------------------------
# Core cutout function
# ---------------------------------------------------------------------------

def make_cutout(imdata, wcs, filtname, src_ID, src_RA, src_Dec, output_dir, size_arcsec=2.0):
    """Create and save a 2-D FITS cutout for a single source.

    Parameters
    ----------
    imdata : numpy.ndarray
        2-D science image array.
    wcs : astropy.wcs.WCS
        WCS object corresponding to *imdata*.
    filtname : str
        Filter name used for directory and filename labelling.
    src_ID : str or int
        Source identifier.
    src_RA : float
        Source right ascension in degrees.
    src_Dec : float
        Source declination in degrees.
    output_dir : str or Path
        Base output directory.  A sub-directory ``{filtname}_cutouts`` is
        created inside it automatically.
    size_arcsec : float, optional
        Side length of the square cutout in arcseconds.  Default is 2.0.

    Returns
    -------
    Path or None
        Path to the saved cutout FITS file, or ``None`` if the cutout was
        skipped (already exists or source is outside the image footprint).
    """
    filt_dir = Path(output_dir) / f"{filtname}_cutouts"
    filt_dir.mkdir(parents=True, exist_ok=True)

    out_image = filt_dir / f"{src_ID}_cutout_{filtname}.fits"

    if out_image.exists():
        logger.debug("Cutout already exists for source %s in %s — skipping.", src_ID, filtname)
        return None

    logger.info("  Creating cutout: source %-12s  filter %s", src_ID, filtname)

    size_deg = size_arcsec / 3600.0
    position = SkyCoord(ra=src_RA, dec=src_Dec, unit="deg")
    size = size_deg * u.deg

    try:
        cutout = Cutout2D(imdata, position, size, wcs=wcs)
        hdu = fits.PrimaryHDU(data=cutout.data, header=cutout.wcs.to_header())
        hdu.writeto(out_image, overwrite=True)
        return out_image
    except ValueError as exc:
        logger.warning("  Skipping source %s (%s): %s", src_ID, filtname, exc)
        return None


# ---------------------------------------------------------------------------
# Field-level processing
# ---------------------------------------------------------------------------

def process_field(src_IDs, src_RAs, src_Decs, im_dir, output_dir,
                  sci_ext="SCI", filter_field_index=1, size_arcsec=2.0):
    """Process all FITS images in a directory and create cutouts for all sources.

    Parameters
    ----------
    src_IDs : array-like
        Source identifiers.
    src_RAs : array-like
        Source right ascensions in degrees.
    src_Decs : array-like
        Source declinations in degrees.
    im_dir : str or Path
        Directory containing the FITS images.
    output_dir : str or Path
        Base output directory for the cutouts.
    sci_ext : str, optional
        HDU extension name for the science data (e.g. ``"SCI"``).
        Default is ``"SCI"``.
    filter_field_index : int, optional
        Index in the underscore-split image filename that gives the filter
        name.  For JADES filenames like ``field_FILTER_...fits`` this is 1.
        Default is 1.
    size_arcsec : float, optional
        Cutout side length in arcseconds.  Default is 2.0.
    """
    im_list = sorted(Path(im_dir).glob("*.fits"))
    if not im_list:
        logger.warning("No FITS files found in %s", im_dir)
        return

    for im_path in im_list:
        parts = im_path.stem.split("_")
        if len(parts) <= filter_field_index:
            logger.warning(
                "Cannot extract filter name from '%s' at index %d — skipping.",
                im_path.name,
                filter_field_index,
            )
            continue

        filtname = parts[filter_field_index]
        logger.info("Processing image: %s  (filter: %s)", im_path.name, filtname)

        with fits.open(im_path, memmap=True) as hdul:
            imdata = hdul[sci_ext].data
            wcs = WCS(hdul[sci_ext].header)

            for src_id, src_ra, src_dec in zip(src_IDs, src_RAs, src_Decs):
                make_cutout(
                    imdata=imdata,
                    wcs=wcs,
                    filtname=filtname,
                    src_ID=src_id,
                    src_RA=src_ra,
                    src_Dec=src_dec,
                    output_dir=output_dir,
                    size_arcsec=size_arcsec,
                )


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="JWST NIRCam cutout creator for JADES GOODS-S and GOODS-N fields.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "catalogue",
        help="Path to the source catalogue CSV file.",
    )
    parser.add_argument(
        "--goods-s-dir",
        required=True,
        help="Directory containing GOODS-S NIRCam FITS images.",
    )
    parser.add_argument(
        "--goods-n-dir",
        required=True,
        help="Directory containing GOODS-N NIRCam FITS images.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Base output directory for cutouts. Overrides config file value.",
    )
    parser.add_argument(
        "--size-arcsec",
        type=float,
        default=None,
        help="Cutout side length in arcseconds. Overrides config file value.",
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Path to a YAML configuration file (see config/example_config.yaml).",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable DEBUG-level logging.",
    )
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    cfg = load_config(args.config)

    # CLI arguments take precedence over config-file values
    size_arcsec = args.size_arcsec if args.size_arcsec is not None else cfg["cutout"]["size_arcsec"]
    output_dir = args.output_dir if args.output_dir is not None else cfg["output"]["dir"]
    sci_ext = cfg["images"]["sci_extension"]
    filter_field_index = cfg["images"]["filter_field_index"]
    cat_cfg = cfg["catalogue"]

    # ------------------------------------------------------------------
    # Load catalogue and split into GOODS-S / GOODS-N subsets
    # ------------------------------------------------------------------
    logger.info("Loading catalogue: %s", args.catalogue)
    src_list = pd.read_csv(args.catalogue)

    gs_mask = src_list[cat_cfg["survey_column"]].isin(cat_cfg["goods_s_surveys"])
    gn_mask = src_list[cat_cfg["survey_column"]].isin(cat_cfg["goods_n_surveys"])
    src_list_gs = src_list[gs_mask]
    src_list_gn = src_list[gn_mask]

    logger.info("GOODS-S sources: %d", len(src_list_gs))
    logger.info("GOODS-N sources: %d", len(src_list_gn))

    id_col = cat_cfg["id_column"]
    ra_col = cat_cfg["ra_column"]
    dec_col = cat_cfg["dec_column"]

    # ------------------------------------------------------------------
    # Process GOODS-S
    # ------------------------------------------------------------------
    if len(src_list_gs) > 0:
        logger.info("=== Processing GOODS-S ===")
        process_field(
            src_IDs=src_list_gs[id_col].to_numpy(),
            src_RAs=src_list_gs[ra_col].to_numpy(),
            src_Decs=src_list_gs[dec_col].to_numpy(),
            im_dir=args.goods_s_dir,
            output_dir=output_dir,
            sci_ext=sci_ext,
            filter_field_index=filter_field_index,
            size_arcsec=size_arcsec,
        )

    # ------------------------------------------------------------------
    # Process GOODS-N
    # ------------------------------------------------------------------
    if len(src_list_gn) > 0:
        logger.info("=== Processing GOODS-N ===")
        process_field(
            src_IDs=src_list_gn[id_col].to_numpy(),
            src_RAs=src_list_gn[ra_col].to_numpy(),
            src_Decs=src_list_gn[dec_col].to_numpy(),
            im_dir=args.goods_n_dir,
            output_dir=output_dir,
            sci_ext=sci_ext,
            filter_field_index=filter_field_index,
            size_arcsec=size_arcsec,
        )

    logger.info("All done!")


if __name__ == "__main__":
    main()
