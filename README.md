# ZVAR-Forced-Phot

Code base for the Caltech Stellar Group's batch forced photometry pipeline.
Original code base can be found here: [przemekmroz/ztf_forced_photometry](https://github.com/przemekmroz/ztf_forced_photometry),
while this forked repo contains the most recent updates, and the pipeline
version that has been used in the most recent forced photometry run
completed in Fall of 2023.

---

## Table of Contents

1. [Overview](#overview)
2. [Repository Structure](#repository-structure)
3. [Pipeline Description](#pipeline-description)
   - [Step 0 – Setup and Authentication](#step-0--setup-and-authentication)
   - [Step 1 – Reference Image Processing](#step-1--reference-image-processing)
   - [Step 2 – Difference Image Photometry](#step-2--difference-image-photometry)
   - [Step 3 – Assembling Epoch Metadata](#step-3--assembling-epoch-metadata)
   - [Step 4 – Building the Matchfile (HDF5 output)](#step-4--building-the-matchfile-hdf5-output)
   - [Step 5 – Zero-Point Calibration (optional)](#step-5--zero-point-calibration-optional)
   - [Step 6 – Light-Curve Metrics (optional)](#step-6--light-curve-metrics-optional)
4. [HDF5 Matchfile Structure](#hdf5-matchfile-structure)
5. [Module Reference](#module-reference)
6. [Utility Programs](#utility-programs)
7. [Running the Pipeline](#running-the-pipeline)
8. [Logging and Error Recovery](#logging-and-error-recovery)
9. [Dependencies](#dependencies)

---

## Overview

ZVAR-Forced-Phot performs **batch PSF forced photometry** on ZTF (Zwicky
Transient Facility) science-minus-reference (difference) images for every
star detected in each ZTF quadrant by the Pan-STARRS1 (PS1) catalog.
Rather than detecting sources anew in every epoch, the pipeline fixes
stellar positions using a reference image and forces a measurement of the
PSF flux at those positions in every difference image.  The resulting
time-series photometry is stored in HDF5 "matchfiles" that encode the
flux, flux error, and quality flag for every star across every observed
epoch.

The pipeline covers two ZTF survey footprints:

| Grid | Field numbers | Typical run script |
|------|---------------|--------------------|
| Primary   | 245 – 881  | `run_phot/run_phot_primary.sh`   |
| Secondary | 1241 – 1897 | `run_phot/run_phot_secondary.sh` |

Each field is subdivided into **16 CCDs × 4 quadrants = 64 quadrant
tiles**, and each is processed independently in all three ZTF filters
(**g**, **r**, **i**).

---

## Repository Structure

```
ZVAR-Forced-Phot/
│
├── psf-phot-pm.py               # Main pipeline driver script
├── psf_dophot.py                # Reference-image photometry via DoPhot
├── build_matchfile.py           # Epoch-file writer and BJD converter
├── metric_suite.py              # Light-curve statistics & ZP calibration
│
├── psf/                         # C-extension: PSF photometry engine
│   ├── psf.c                    #   Core PSF flux estimator (CFITSIO)
│   ├── psf_wrapper.c            #   Python/C-API wrapper
│   ├── setup.py                 #   Build script (distutils)
│   ├── build_module.sh          #   Convenience build script
│   └── PSF_Phot_PM.so           #   Pre-built shared library
│
├── matchfiles/                  # Matchfile tools (C and Python)
│   ├── build_matchfile.c        #   Writes HDF5 matchfiles
│   ├── matchfile_pm.c / .h      #   HDF5 read library
│   ├── getobj.c / getobj.py     #   Extract individual light curves
│   ├── AOV_matchfile.cu         #   GPU-based Analysis of Variance (CUDA)
│   ├── AOV_matchfile_single_prec.cu  # Single-precision CUDA variant
│   ├── get_scatter.c            #   Per-star photometric scatter
│   ├── read_data.c              #   Low-level HDF5 reader
│   ├── radec2field.py           #   RA/Dec → field/CCD/quad resolver
│   ├── utils.py                 #   Geometric helper functions
│   └── ZTF_Fields.txt           #   ZTF field center coordinates
│
├── dophot/                      # DoPhot PSF fitting binary
│   ├── dophot                   #   Compiled DoPhot v4 binary
│   └── param_default_c          #   Default DoPhot parameter template
│
├── cals/                        # Zero-point calibration
│   ├── calibrate_to_ZTF_system.py         # ZP correction (per-quadrant)
│   ├── calibrate_to_ZTF_system_matchfile.py  # ZP correction (from matchfile)
│   ├── get_zp_corr.sh           #   Batch ZP-correction runner
│   ├── zero_point_corrections_primary.txt   # Saved g/r/i ZP offsets (primary)
│   └── zero_point_corrections_secondary.txt # Saved g/r/i ZP offsets (secondary)
│
├── run_phot/                    # Batch-run shell scripts
│   ├── run_phot_primary.sh      #   Primary-grid run
│   ├── run_phot_secondary.sh    #   Secondary-grid run
│   └── run_phot_rerun.sh        #   Re-run failed quadrants
│
├── logs/                        # Logging utilities
│   ├── checklog.py              #   Parse logs, extract failed quadrants
│   ├── checklog.sh              #   Batch log checker
│   └── log_files/               #   Example log files per field
│
├── sBLS/                        # Sparse Box Least Squares
│   ├── sbls.c / sbls.h          #   BLS implementation
│   └── period_find.c            #   Period-search driver
│
├── syncfiles/                   # Data synchronization
│   ├── copy_matchfiles_to_schoty.sh  # rsync matchfiles to schoty server
│   └── check_rsynclog.py        #   Verify rsync completed successfully
│
└── docs/                        # Documentation (LaTeX + PDF)
    ├── ZTF_forcedphot_docs.*    #   Pipeline technical description
    ├── ZTF_forcedphot_runpipe.* #   How to run the pipeline
    └── ZTF_forcedphot_updates.* #   Changelog / update notes
```

---

## Pipeline Description

The main entry point is `psf-phot-pm.py`.  It is invoked with the
field number, CCD number, quadrant number, and filter, plus optional
Julian Date limits:

```bash
python3 psf-phot-pm.py <field> <ccd> <quad> <filter> \
    [--jd_start JD_START] [--jd_end JD_END] [--logfile LOGFILE]

# Example
python3 psf-phot-pm.py 760 14 3 r --jd_start 2458200 --jd_end 2458600
```

The diagram below summarises the data flow through the pipeline:

```
ZTF Archive (IRSA)           Kowalski / Gloria / Melman databases
        │                              │
        │  reference image             │  PS1 DR1 star list
        ▼                              ▼
┌─────────────────────────────────────────────────┐
│  Step 1 – Reference Image Processing            │
│  psf_dophot.py → DoPhot binary                  │
│  Outputs: starlist.txt, dophot_output.txt,       │
│           params.json, ps_sources.csv            │
└──────────────────────┬──────────────────────────┘
                       │  fixed source positions
        ┌──────────────▼──────────────────────────┐
        │  ZTF Archive (IRSA)                      │
        │  difference images + PSF + science imgs  │
        └──────────────┬──────────────────────────┘
                       │
┌──────────────────────▼──────────────────────────┐
│  Step 2 – Difference Image Photometry           │
│  psf-phot-pm.py → PSF_Phot_PM.so (C extension) │
│  Outputs: *_psf.cat binary files (per epoch)    │
└──────────────────────┬──────────────────────────┘
                       │
┌──────────────────────▼──────────────────────────┐
│  Step 3 – Epoch Metadata Assembly               │
│  build_matchfile.py  →  epochs.txt              │
│  (JD → BJD conversion; image quality flags)     │
└──────────────────────┬──────────────────────────┘
                       │
┌──────────────────────▼──────────────────────────┐
│  Step 4 – Build HDF5 Matchfile                  │
│  matchfiles/build_matchfile (compiled C)         │
│  Output: data_FFFF_CC_Q_ZZ.h5                   │
└──────────────────────┬──────────────────────────┘
                       │
           ┌───────────┴───────────┐
           ▼                       ▼
┌──────────────────┐   ┌──────────────────────────┐
│  Step 5 (opt.)   │   │  Step 6 (opt.)            │
│  ZP Calibration  │   │  Light-Curve Metrics      │
│  cals/           │   │  metric_suite.py           │
└──────────────────┘   └──────────────────────────┘
```

### Step 0 – Setup and Authentication

Before doing any science, the script establishes connections to three
Kowalski-family databases using credentials read from `/data/passwd.json`:

| Instance | Hostname | Purpose |
|----------|----------|---------|
| **Kowalski** | `kowalski.caltech.edu` | ZTF operations database (`ZTF_ops`), Gaia EDR3, PS1 DR1 |
| **Gloria**   | `gloria.caltech.edu`   | ZTF source-feature catalog (`ZTF_source_features_DR5`) for SCoPE metrics |
| **Melman**   | `melman.caltech.edu`   | ZTF exposure metadata (`ZTF_exposures_*`) |

The pipeline also uses IRSA (NASA/IPAC Infrared Science Archive) to
download images, authenticating via a cookie stored at `~/.ztf_cookie`.

---

### Step 1 – Reference Image Processing

**Relevant code:** `psf_dophot.py`, `psf-phot-pm.py`
(`check_if_reference_images_are_processed`, `download_reference_images`,
`process_ref_images`)

#### 1a. Check for existing work

The function `check_if_reference_images_are_processed` looks for three
files in the results directory:

| File | Contents |
|------|----------|
| `*_starlist.txt` | PS1 source ID, RA, Dec, x, y (pixel) |
| `*_dophot_input.txt` | DoPhot input (position + Gaussian shape) |
| `params.json` | Image metadata: ZP, FWHM, gain, sky, saturation, WCS |

If all three exist, this step is skipped entirely.

#### 1b. Download the reference image

A single co-added reference image per quadrant/filter is downloaded from
IRSA at:

```
https://irsa.ipac.caltech.edu/ibe/data/ztf/products/ref/
  {000|001}/field{FFFFFF}/z{g|r|i}/ccd{CC}/q{Q}/
  ztf_{FFFFFF}_z{g|r|i}_c{CC}_q{Q}_refimg.fits
```

Primary fields (< 1000) use prefix `000`; secondary fields use `001`.

#### 1c. Query the PS1 catalog

`get_PS1_stars` queries Kowalski's `PS1_DR1` catalog within a 0.64°
cone of the image centre, selecting only sources with ≥ 5 detections.
The resulting PS1 source IDs, RA, Dec, g-band and r-band PSF magnitudes
are written to `ps_sources.csv` and `*_starlist.txt`.

#### 1d. Run DoPhot

DoPhot v4 (`dophot/dophot`) is invoked with a per-quadrant parameter
file that sets:

- `PSFTYPE = PGAUSS` – position-fixed Gaussian PSF model
- `FIXPOS = YES` – centroids are held at the PS1 positions
- `FWHM`, `SKY`, `EPERDN`, `RDNOISE` – from the reference image header
- Detection thresholds (`THRESHMIN`, `THRESHMAX`), saturation
  (`CTPERSAT`, `ITOP`), and various I/O file paths

Output is written to `*_dophot_output.txt`.  This file records, for
every input PS1 star, the fitted magnitude, magnitude uncertainty, and
object type (1 = point source).

---

### Step 2 – Difference Image Photometry

**Relevant code:** `psf-phot-pm.py` (`download_ZTF_images`,
`process_single_image_psf`, `process_images_psf`), `psf/psf.c`,
`psf/psf_wrapper.c`

#### 2a. Query the ZTF operations database

`get_image_list` queries `ZTF_ops` on Kowalski for all exposures of the
requested field/CCD/quad within the requested Julian Date range.  The
following images are filtered out:

- Exposures with `pid = 0` (non-science)
- Public-program exposures (`pid = 1`) obtained after the public-data
  cutoff date (`PUBLIC_DATE = 2460071.5` ≈ 2023-05-07, DR18)

#### 2b. Download images

For each qualifying frame, three FITS files are downloaded from IRSA in
parallel (default 16 threads):

| Image type | IRSA suffix |
|------------|-------------|
| Difference image (gzip-compressed) | `_scimrefdiffimg.fits.fz` |
| Science image | `_sciimg.fits` |
| PSF image (DAO) | `_sciimgdaopsfcent.fits` |

Already-downloaded files are skipped.

#### 2c. Quality filtering

Before performing photometry, `process_single_image_psf` reads the
difference-image FITS header and applies per-filter quality cuts:

| Parameter | ZTF_g | ZTF_r | ZTF_i |
|-----------|-------|-------|-------|
| `INFOBITS` | 0 | 0 | 0 |
| `STATUS` | ≥ 1 | ≥ 1 | ≥ 1 |
| Seeing (arcsec) | < 3.0 | < 3.0 | < 3.0 |
| Exposure time (s) | ≤ 30 | ≤ 30 | ≤ 30 |
| `NMATCHES` (PS1 cal stars) | ≥ 80 | ≥ 120 | ≥ 100 |
| ZP scatter (`MAGZPRMS`) | ≤ 0.06 | ≤ 0.05 | ≤ 0.06 |
| ZP range | 25.80 – 26.70 (airmass-corrected) | 25.80 – 26.65 | 25.80 – 26.00 |

Images failing any cut are skipped without writing an output file.

#### 2d. PSF flux estimation

The C extension `PSF_Phot_PM` (`psf/psf.c`) performs the actual
measurement.  For each star at pixel position (x, y) the algorithm:

1. Reads a `(2·stamp+1)²` pixel postage stamp (default `stamp=3`, so
   7 × 7 pixels) from the difference image.
2. Looks up the PSF model value at each pixel by bilinear interpolation
   of the PSF image centred at the star's sub-pixel position.
3. Computes the optimal (inverse-variance weighted) flux estimate:

   ```
   flux  = Σ [ T_ki · D_ki / σ²_ki ] / Σ [ T²_ki / σ²_ki ]
   ferr  = 1 / sqrt( Σ [ T²_ki / σ²_ki ] )
   ```

   where `D_ki` is the difference-image pixel value, `T_ki` is the PSF
   template value, and `σ²_ki = science_pixel / gain` is the pixel
   variance estimated from the science image.

4. Returns `NaN` if fewer than half the stamp pixels are valid, or if
   the denominator is non-positive.

The flux is then rescaled to a **fiducial zero-point** of 26.275 mag:

```
flux_fid = flux · 10^{ 0.4 · (ZP_fid − ZP_image) }
```

Results for all stars in one epoch are written to a binary file
`*_psf.cat` (two arrays of 32-bit floats: `flux[N_stars]` then
`ferr[N_stars]`).

---

### Step 3 – Assembling Epoch Metadata

**Relevant code:** `psf-phot-pm.py` (`process_images_psf`),
`build_matchfile.py` (`save_to_file`, `JD_to_BJD`)

After photometry, each successfully processed epoch contributes one row
to the `epochs.txt` file with 21 space-separated columns:

| Column | Name | Description |
|--------|------|-------------|
| 1 | `jd` | Mid-exposure Julian Date |
| 2 | `bjd` | Barycentric Julian Date (TDB), corrected for light travel time to Palomar (lat 33°21′26″, lon −116°51′36″, elev 1712 m) |
| 3 | `fid` | Filter ID (1=g, 2=r, 3=i) |
| 4 | `exptime` | Exposure time (s) |
| 5 | `pid` | Program ID (1=public, 2=private) |
| 6 | `field` | ZTF field number |
| 7 | `ccd` | CCD number (1–16) |
| 8 | `quad` | Quadrant number (1–4) |
| 9 | `imstat` | Image status flag |
| 10 | `infobits` | Quality infobits |
| 11 | `seeing` | Seeing FWHM (arcsec) |
| 12 | `magzp` | Magnitude zero-point |
| 13 | `magzpunc` | ZP uncertainty |
| 14 | `magzprms` | ZP RMS scatter |
| 15 | `clrcoeff` | Color coefficient |
| 16 | `clrcounc` | Color coefficient uncertainty |
| 17 | `limmag` | 5σ limiting magnitude |
| 18 | `airmass` | Airmass |
| 19 | `nps1matches` | Number of PS1 calibration stars |
| 20 | `qcomment` | Quadrant comment string |
| 21 | `result_file` | Path to the binary `*_psf.cat` file |

---

### Step 4 – Building the Matchfile (HDF5 output)

**Relevant code:** `matchfiles/build_matchfile.c`,
`psf-phot-pm.py` (`process_quad` → subprocess call to `build_matchfile`)

The compiled C program `matchfiles/build_matchfile` reads the text
files produced in Steps 1–3 and writes a single compressed HDF5 file
per quadrant/filter combination.  The file is named:

```
data_{FFFF}_{CC}_{Q}_{ZZ}.h5
```

where `FFFF` = zero-padded field, `CC` = CCD, `Q` = quadrant,
`ZZ` = filter code (`zg`, `zr`, or `zi`).

See the [HDF5 Matchfile Structure](#hdf5-matchfile-structure) section
below for the full schema.  Datasets are stored with gzip compression
(level 1) and a chunk size of 1024 elements.

---

### Step 5 – Zero-Point Calibration (optional)

**Relevant code:** `cals/calibrate_to_ZTF_system.py`,
`cals/get_zp_corr.sh`,
`cals/zero_point_corrections_primary.txt`,
`cals/zero_point_corrections_secondary.txt`,
`metric_suite.py` (`mag_corr`)

The DoPhot magnitudes in the reference image use that image's own
zero-point.  A systematic offset relative to the ZTF photometric system
is calculated by cross-matching DoPhot-measured stars against median
magnitudes from the ZTF source catalog (`ZTF_sources_20230821` on
Gloria).  The offset is the median of `(ZTF_mag − DoPhot_mag)` over all
matched point sources within 1 arcsec.

The batch script `cals/get_zp_corr.sh` iterates over all
field/CCD/quadrant combinations and stores results as:

```
FFFF CC Q  ZP_g  N_g  ZP_r  N_r  ZP_i  N_i
```

At runtime, `metric_suite.mag_corr` applies these corrections when
computing magnitudes from the forced-photometry fluxes.

---

### Step 6 – Light-Curve Metrics (optional)

**Relevant code:** `metric_suite.py` (`calc_lc_metrics`,
`get_scope_features`)

For a configurable set of ZTF fields (`METRIC_FIELDS` in
`psf-phot-pm.py`; empty by default), the pipeline computes a suite of
variability statistics over each star's light curve.

#### 6a. SCoPE cross-match

`get_scope_features` queries the `ZTF_source_features_DR5` catalog on
Gloria (paginating in batches of up to 10 000 entries), then
cross-matches within 1 arcsec to the quadrant's starlist.  For each
matched star the following quantities are retrieved:

- Median magnitude
- Best-fit period (from the SCoPE survey)
- Period significance

#### 6b. Variability statistics

`calc_lc_metrics` loads all `*_psf.cat` binary files, converts
differential fluxes to magnitudes using `flux2mag` (adding the DoPhot
reference-image magnitude, propagating uncertainties, and applying ZP
corrections), then applies quality cuts (infobits = 0, seeing < 3″,
ZP scatter < 0.06, airmass < 2.5, ≥ 80 calibration stars), and
computes 12 statistics for each star with at least 20 valid points:

| Metric | Description |
|--------|-------------|
| `stdev` | Standard deviation of magnitudes |
| `pp_scat` | RMS point-to-point scatter |
| `ms_successive_diff` | Mean square successive difference |
| `skew` | Normalized skewness (error-weighted) |
| `skewP` | Pearson skewness |
| `kurtosis` | Normalized kurtosis |
| `norm_ev` | Normalized excess variance |
| `von_neumann` | Von Neumann ratio (MSSD / variance) |
| `interq_range` | Inter-quartile range |
| `welch_i` | Welch/Stetson I variability index |
| `stetson_J` | Stetson J index |
| `stetson_K` | Stetson K index |

Results are written to `lc_metrics.txt` for optional later ingestion
into the matchfile.

---

## HDF5 Matchfile Structure

Each matchfile contains a single top-level group `/data` with three
datasets and one attribute:

```
/data  (group)
│   └── attribute: magzp  (float32) – reference-image ZP magnitude
│
├── /data/sources  (compound dataset, N_stars rows)
│       gaia_id      uint64  – PS1 source ID (stored in Gaia_ID field)
│       ra           float32 – Right ascension (deg)
│       decl         float32 – Declination (deg)
│       mag_ref      float32 – Reference-image magnitude (DoPhot)
│       mag_err_ref  float32 – Reference-image magnitude error
│       objtype      uint16  – DoPhot object type (1 = point source)
│
├── /data/exposures  (compound dataset, N_exp rows)
│       jd           float64 – Julian Date (mid-exposure)
│       bjd          float64 – Barycentric Julian Date (TDB)
│       filterid     uint16  – Filter (1=g, 2=r, 3=i)
│       exptime      float32 – Exposure time (s)
│       pid          uint16  – Program ID
│       field        uint16  – ZTF field number
│       ccd          uint16  – CCD number (1–16)
│       quad         uint16  – Quadrant number (1–4)
│       imstat       uint16  – Image status
│       infobits     uint64  – Quality bitmask
│       seeing       float32 – Seeing FWHM (arcsec)
│       mzpsci       float64 – Science-image ZP
│       mzpsciunc    float64 – ZP uncertainty
│       mzpscirms    float64 – ZP RMS scatter
│       clrco        float64 – Color coefficient
│       clrcounc     float64 – Color coefficient uncertainty
│       maglim       float32 – 5σ limiting magnitude
│       airmass      float32 – Airmass
│       nps1matches  uint64  – Number of PS1 calibration stars
│
└── /data/sourcedata  (compound dataset, N_stars × N_exp rows)
        flux         float32 – Forced-PSF flux (fiducial ZP 26.275)
        flux_err     float32 – Flux uncertainty
        flag         uint16  – Quality flag (0=good, 1=NaN measurement)
```

The `sourcedata` array is stored in row-major order: element
`[star_index * N_exp + epoch_index]` gives the measurement of star
`star_index` at epoch `epoch_index`.

---

## Module Reference

### `psf-phot-pm.py` — Main pipeline driver

| Function | Description |
|----------|-------------|
| `get_image_list` | Query `ZTF_ops` on Kowalski for exposures in a field/CCD/quad/JD range |
| `get_cookie` | Fetch and cache an IRSA authentication cookie |
| `load_file` | Download a single file from IRSA using streaming HTTP |
| `download_data` | Parallel download of a list of URLs (multiprocessing Pool) |
| `download_ZTF_images` | Download difference, PSF, and science images for a quadrant |
| `get_image_info` | Extract 17 quality and calibration parameters from a FITS header |
| `process_single_image_psf` | Quality-filter one epoch and run forced PSF photometry |
| `process_images_psf` | Run `process_single_image_psf` in parallel over many images |
| `get_image_info_preprocessed` | Reconstruct epoch metadata from Kowalski/IRSA when the original header is unavailable |
| `download_reference_images` | Download a quadrant reference image from IRSA |
| `process_ref_images` | Run DoPhot on a reference image (via `psf_dophot.get_phot_ref_image`) |
| `check_if_reference_images_are_processed` | Check for the three sentinel files produced by Step 1 |
| `check_if_diffim_photometry_exists` | Inventory existing `*_psf.cat` files and reconcile with `epochs.txt` |
| `process_quad` | Top-level orchestrator: runs Steps 1–6 for one field/CCD/quad/filter |

### `psf_dophot.py` — Reference-image photometry

| Function | Description |
|----------|-------------|
| `get_info_from_ref_image` | Read header metadata (ZP, FWHM, WCS, saturation) from a reference FITS file |
| `radec_to_xy` | Convert RA/Dec to pixel coordinates using Astropy WCS |
| `get_Gaia_stars` | Cone search of Gaia EDR3 on Kowalski |
| `get_PS1_stars` | Cone search of PS1 DR1 on Kowalski (filters by `nDetections ≥ 5`) |
| `create_dophot_parameter_input_file` | Write a per-quadrant DoPhot parameter file |
| `get_phot_ref_image` | End-to-end reference-image pipeline: PS1 query → starlist → DoPhot → `params.json` |

### `build_matchfile.py` — Epoch file and BJD conversion

| Function | Description |
|----------|-------------|
| `JD_to_BJD` | Convert JD array to BJD(TDB) accounting for light travel time to Palomar Observatory |
| `save_to_file` | Sort epochs by JD, add BJD column, write `epochs.txt`, and return the HDF5 output filename |

### `metric_suite.py` — Light-curve statistics and calibration

| Function | Description |
|----------|-------------|
| `stdev`, `pp_scat`, `ms_successive_diff`, `skew`, `skewP`, `kurtosis`, `norm_ev`, `von_neumann`, `interq_range`, `welch_i`, `stetson_J`, `stetson_K` | Per-star variability statistics (vectorised over all stars simultaneously) |
| `mag_corr` | Apply a stored zero-point correction from `cals/zero_point_corrections_*.txt` |
| `flux2mag` | Convert forced-photometry differential flux to total magnitude (difference flux + reference flux from DoPhot) and apply ZP correction |
| `get_scope_features` | Query Gloria for SCoPE periodicity features and cross-match to the starlist |
| `calc_lc_metrics` | Load all photometry, convert to magnitudes, apply quality cuts, compute all 12 statistics, write `lc_metrics.txt` |

### `matchfiles/utils.py` — Geometric helpers

| Function | Description |
|----------|-------------|
| `ang_dist` | Great-circle angular distance (radians) |
| `ortographic_projection` | Project RA/Dec onto the tangent plane of a field centre |
| `inside_polygon` | Determine which ZTF CCD and quadrant contains a given tangent-plane coordinate |
| `read_fields` | Load the ZTF field grid from `ZTF_Fields.txt` |
| `get_field_id` | Find all field/CCD/quad tiles that contain a given RA/Dec |
| `get_PS1_stars` | Cone search of PS1 for nearby sources (used by `radec2field.py`) |

---

## Utility Programs

### `matchfiles/getobj.py` (and compiled `matchfiles/getobj`)

Extract the light curve of a single star from a matchfile:

```bash
python3 matchfiles/getobj.py  DATA_FILE  OBJECT_ID

# Output columns: BJD  flux  flux_err  flag  pid
```

The compiled C version `matchfiles/getobj` provides the same
functionality at higher speed.

### `matchfiles/radec2field.py`

Look up which field/CCD/quadrant and PS1 source ID correspond to a
given sky position:

```bash
python3 matchfiles/radec2field.py  RA_DEG  DEC_DEG  [--radius ARCSEC]
```

### `matchfiles/AOV_matchfile.cu`

GPU-accelerated Analysis of Variance (AOV) period search over an entire
matchfile.  Compiled with:

```bash
nvcc -arch=sm_60 -o AOV AOV_matchfile.cu matchfile_pm.c \
    -I /path/to/hdf5/include -L /path/to/hdf5/lib -lhdf5
```

Usage: `./AOV  MATCHFILE  PERIOD_MIN  PERIOD_MAX`

### `sBLS/` — Sparse BLS period finder

A pure-C implementation of the sparse Box Least Squares algorithm for
detecting transit-like signals.  Build with `make` and run
`period_find`.

---

## Running the Pipeline

### Prerequisites

- Python ≥ 3.8 with packages: `numpy`, `pandas`, `astropy`,
  `penquins`, `requests`, `h5py`, `scipy`, `matplotlib`
- CFITSIO (for compiling the PSF C extension)
- HDF5 C library (for `matchfiles/build_matchfile`)
- Access credentials in `/data/passwd.json`:

  ```json
  {
    "kowalski":  {"username": "...", "password": "..."},
    "gloria":    {"token": "..."},
    "irsa":      {"username": "...", "password": "..."}
  }
  ```

### Build the PSF module

```bash
cd psf
python3 setup.py build_ext --inplace
# Or use the convenience script:
bash build_module.sh
```

### Build the matchfile tools

```bash
cd matchfiles
gcc -o build_matchfile build_matchfile.c matchfile_pm.c \
    -I /usr/local/hdf5/include -L /usr/local/hdf5/lib -lhdf5 -lm
```

### Run a single quadrant

```bash
python3 psf-phot-pm.py 760 14 3 r --jd_start 2458200 --jd_end 2460157
```

### Run a full field range (batch)

```bash
# Primary fields (245–881):
bash run_phot/run_phot_primary.sh

# Secondary fields (1241–1897):
bash run_phot/run_phot_secondary.sh
```

Each quadrant processes g, r, and i filters sequentially, cleans up
large temporary FITS files between filters, and logs progress to
`logs/log_files/pipeline_{FFFF}.log`.

---

## Logging and Error Recovery

Every invocation of `psf-phot-pm.py` writes structured log entries to
a per-field log file.  The format of a successful entry is:

```
TIMESTAMP:INFO:f{FFFFFF}_c{CC}_q{Q}_z{f} SUCCESS
   # Downloading data      : X.XX sec
   # Doing photometry      : X.XX sec
   ...
   # Number of stars: NNNN
   # Number of images processed: MMMM
```

Failed quadrants are tagged `FAILED` or `MetricError`.

**Identify and re-run failures:**

```bash
# 1. Parse all log files for the processed field range
bash logs/checklog.sh
# → generates logs/failed_quads.txt

# 2. Re-run only the failed quadrants
bash run_phot/run_phot_rerun.sh
```

---

## Dependencies

| Package / Tool | Purpose |
|----------------|---------|
| `numpy`, `pandas` | Array and DataFrame operations |
| `astropy` | FITS I/O, WCS, time conversions |
| `penquins` | Kowalski/Gloria/Melman database client |
| `requests` | HTTP downloads from IRSA |
| `h5py` | Reading HDF5 matchfiles in Python |
| `scipy` | `scipy.stats.iqr` for the inter-quartile range metric |
| `matplotlib` | (Imported in metric_suite; used for optional diagnostic plots) |
| CFITSIO | C library used by `psf/psf.c` to read FITS images |
| HDF5 C library | Used by `matchfiles/build_matchfile.c` to write matchfiles |
| DoPhot v4 | Compiled binary in `dophot/dophot` for reference-image PSF photometry |
| CUDA toolkit | Required only for `matchfiles/AOV_matchfile.cu` (GPU period search) |
