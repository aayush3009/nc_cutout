# nc_cutout

**JWST NIRCam postage-stamp cutout creator for JADES GOODS-S and GOODS-N data**

`nc_cutout` takes a source catalogue (CSV) and a directory of NIRCam FITS images and produces small, WCS-correct postage-stamp cutouts for every source in every available filter.

---

## Features

- Processes both **GOODS-S** and **GOODS-N** NIRCam mosaics in a single run
- Skips cutouts that already exist so interrupted runs can be safely resumed
- Configurable via **command-line arguments** and/or a **YAML config file**  
  (CLI arguments always win)
- Memory-efficient: reads FITS files with `memmap=True` and only extracts pixels on demand via `astropy.nddata.Cutout2D`
- Structured logging with an optional `--verbose` flag for debug output

---

## Requirements

- Python ≥ 3.9
- [astropy](https://www.astropy.org/) ≥ 5.0
- [pandas](https://pandas.pydata.org/) ≥ 1.4
- [PyYAML](https://pyyaml.org/) ≥ 6.0

---

## Installation

Clone the repository and install in editable mode:

```bash
git clone https://github.com/aayush3009/nc_cutout.git
cd nc_cutout
pip install -e .
```

Or install directly from GitHub:

```bash
pip install git+https://github.com/aayush3009/nc_cutout.git
```

---

## Quick start

```bash
nc_cutout /path/to/catalogue.csv \
    --goods-s-dir /Volumes/data/JADES-NIRCam/GOODS-S-NC-v0.9/ \
    --goods-n-dir /Volumes/data/JADES-NIRCam/GOODS-N-NC-v0.7/ \
    --output-dir ./cutouts \
    --size-arcsec 2.0
```

Cutouts are saved under `./cutouts/{FILTER}_cutouts/{ID}_cutout_{FILTER}.fits`.

---

## Using a config file

Copy and edit the provided example:

```bash
cp config/example_config.yaml my_config.yaml
# edit my_config.yaml to match your column names, survey tags, etc.

nc_cutout catalogue.csv \
    --goods-s-dir /path/to/goods-s/ \
    --goods-n-dir /path/to/goods-n/ \
    --config my_config.yaml
```

Any setting in the config file can still be overridden on the command line.

---

## Catalogue format

The script expects a **CSV file** with at least these columns (names are configurable):

| Column | Default name | Description |
|--------|-------------|-------------|
| Source ID | `ID_1` | Unique identifier — used in output filenames |
| RA | `RA` | Right ascension in **degrees** |
| Dec | `DEC` | Declination in **degrees** |
| Survey | `SURVEY` | Survey name used to route source to GOODS-S or GOODS-N |

---

## Image file naming convention

Filter names are extracted from the FITS filename by splitting on `_` and taking the element at `filter_field_index` (default 1). This matches the standard JADES naming scheme:

```
<field>_<FILTER>_<version>_<...>.fits
          ^
          index 1
```

If your filenames follow a different convention, adjust `filter_field_index` in your config file.

---

## Command-line reference

```
usage: nc_cutout [-h] --goods-s-dir GOODS_S_DIR --goods-n-dir GOODS_N_DIR
                 [--output-dir OUTPUT_DIR] [--size-arcsec SIZE_ARCSEC]
                 [--config CONFIG] [--verbose]
                 catalogue

positional arguments:
  catalogue                     Path to the source catalogue CSV file.

options:
  --goods-s-dir GOODS_S_DIR     Directory containing GOODS-S NIRCam FITS images.
  --goods-n-dir GOODS_N_DIR     Directory containing GOODS-N NIRCam FITS images.
  --output-dir  OUTPUT_DIR      Base output directory (default: ./cutouts).
  --size-arcsec SIZE_ARCSEC     Cutout side length in arcseconds (default: 2.0).
  --config      CONFIG          Path to a YAML configuration file.
  --verbose                     Enable DEBUG-level logging.
```

---

## Python API

The core functions can also be used directly in notebooks or scripts:

```python
from nc_cutout import make_cutout, process_field
from astropy.io import fits
from astropy.wcs import WCS

with fits.open("my_image.fits", memmap=True) as hdul:
    imdata = hdul["SCI"].data
    wcs    = WCS(hdul["SCI"].header)

make_cutout(
    imdata=imdata,
    wcs=wcs,
    filtname="F200W",
    src_ID="GS-12345",
    src_RA=53.1234,
    src_Dec=-27.5678,
    output_dir="./my_cutouts",
    size_arcsec=2.0,
)
```

---

## Contributing

Pull requests and issues are welcome!  
Please open an issue before submitting large changes.

---

## License

MIT — see [LICENSE](LICENSE) for details.
