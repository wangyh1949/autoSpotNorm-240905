# autoSpotNorm

`autoSpotNorm` is a MATLAB workflow for bacterial single-particle-tracking (SPT) experiments. It combines cell meshes produced by [Oufti](https://oufti.org/) with trajectories produced by u-track, assigns tracks to cells, corrects a global image offset, calculates each spot's normalized location within a cell, and derives diffusion and step-based measurements.

The workflow defines a spot position as `[LNorm, xNorm]`:

- `LNorm` is the normalized position along the cell long axis (`0` to `1`).
- `xNorm` is the normalized transverse position relative to the cell centerline (approximately `-1` to `1`).

## Requirements

- MATLAB, with the functions used by this project available (including `boundary` and `exportgraphics`).
- Oufti mesh output containing `cellList` and `cellListN`.
- A u-track result named `Channel_1_tracking_result.mat`, containing `tracksFinal`.
- The repository root and its `Function` directory on the MATLAB path. From MATLAB, run:

```matlab
addpath(genpath('path-to-autoSpotNorm'))
```

## Input layout and naming

Select a folder for one experiment day whose name follows this convention:

```text
YYMMDD-SK### [optional experiment suffix]
```

For example, `230712-SK187 Suc`. The date, strain, and optional suffix are parsed from this folder name and used in output names.

Within that folder, the scripts expect the following:

```text
230712-SK187 Suc/
  BF001m.mat
  BF001.tif
  Tracking001/
    TrackingPackage/tracks/Channel_1_tracking_result.mat
```

Important conventions:

- Mesh files are discovered with `*BF*m.mat` and are expected to use a three-digit movie number, such as `BF001m.mat`.
- The corresponding TIFF image is inferred from the mesh name (`BF001.tif`).
- The tracking folder is located by matching `*rack*` followed by the movie number (for example, `Tracking001`).
- Each mesh file must contain Oufti's `cellList` and `cellListN`; each tracking result must contain u-track's `tracksFinal`.

## Single-day workflow

Use `autoSpotNorm.m` to process all mesh/tracking pairs in a selected experiment folder and combine their results.

Before running it, edit the settings at the top of the script:

- `dateAdd`: label appended to the experiment date in result names.
- `pixelSize`: camera pixel size in metres (default: `160e-9`).
- `timeStep`: frame interval in seconds (default: `21.742e-3`).
- `imgSaveFlag`: whether to save mesh/track overlay PNG files.
- `goodCellFlag`: `true` accepts all Oufti cells as isolated; `false` prompts for non-isolated cells.
- `autoCombineMovieFlag`: `true` combines all eligible cells; `false` prompts for the cells to include per movie.
- `varPath`: root directory for analysis output. Change this from the author-specific default.

Then run:

```matlab
autoSpotNorm
```

The script prompts for a representative mesh file, switches to that folder, processes every matching mesh file, asks whether the strain is cytoplasmic, and uses this to set the mostly-filled-cell threshold used for image-shift correction.

### What happens for each movie

`spotNorm_yh.m` performs the per-movie processing:

1. Loads the Oufti mesh and u-track trajectories.
2. Converts u-track coordinates to `tracksCoordXY` and assigns each track to the cell containing most of its spots.
3. Generates an initial mesh/track overlay.
4. Estimates a global XY shift from cells whose track coverage occupies enough of the cell area, then applies it to all tracks. A movie without at least three suitable cells is marked as `badMovie`.
5. Retains tracks in selected isolated cells and calculates `spotPosNorm` for every in-cell spot.
6. Saves a per-movie `_Variables.mat` file.

## Multi-day workflow

Use `combineTF.m` after individual experiment days have been processed. It finds the relevant per-movie `_Variables.mat` files under `varPath`, lets you select experiment dates, combines their tracks, and then performs the same diffusion, matrix, and step analyses.

Set `varPath` at the top of `combineTF.m`, then run:

```matlab
combineTF
```

Enter the strain identifier(s) and, when offered, the experiment-day indices to combine. Supply the requested combined-date label when prompted.

## Output

All outputs are MATLAB MAT files under the configured `varPath`.

| Output | Produced by | Contents |
| --- | --- | --- |
| `Variables/<strain>/<date>/*_Variables.mat` | `spotNorm_yh.m` | Per-movie meshes, tracks, shift diagnostics, cell-selection flags, and `tracksFinal_InIsolatedCell`. |
| `Track Analysis Archive/Single Day/<name>.mat` | `autoSpotNorm.m` | Combined tracks and cell metadata before diffusion fields are added. |
| `Matrix/Single Day/<name>.mat` | `autoSpotNorm.m` | Plot-ready per-track location, diffusion, MSD, cell-region, and metadata arrays. |
| `Track Analysis Archive/<strain>/<name>.mat` | `combineTF.m` | Combined multi-day tracks and cell metadata before diffusion fields are added. |
| `tracksFinal/<name> All.mat` | `combineTF.m` | Combined multi-day tracks with diffusion-analysis fields. |
| `Matrix/<name>.mat` | `combineTF.m` | Plot-ready multi-day per-track arrays. |
| `Steps/<name>.mat` | `combineTF.m` | Per-step displacement, longitudinal displacement, normalized position, orientation, and radius-of-gyration data. |
| `Mesh & Tracks Image/*.png` | `spotNorm_yh.m` | Initial overlays; shifted overlays are in `auto shiftBack/` when image saving is enabled. |

## Analysis fields and arrays

`tracksFinal` is augmented with:

- `tracksCoordXY`: XY spot coordinates in pixels (shift-corrected when applicable).
- `cellNumber`: original Oufti cell index.
- `ModCellNum`: consecutively renumbered cell index after combining movies.
- `spotPosNorm`: `N x 2` `[LNorm, xNorm]` values for a track's spots.
- `origin`, `badMovie`, and `filled`: source movie and quality/coverage flags.
- `MSD`, `Diff`, `LocErr`, `Alpha`, and `Dalpha`: added by diffusion analysis.

Diffusion analysis uses time-averaged MSD and fits the first three lag times. `Diff` is the apparent diffusion coefficient, `LocErr` is localization error, and `Alpha`/`Dalpha` come from a power-law MSD fit. Ensemble MSD and time-averaged MSD arrays retain up to 50 lag times.

The matrix files contain compact, plot-ready forms of these results, including first/end/min/max `LNorm`, first/end and first-four-frame `xNorm`, diffusion metrics, and masks for tracks or spots outside the pole regions. Pole regions are defined as one-half cell width at each end of the long axis.

## Supporting scripts

| Script | Purpose |
| --- | --- |
| `spotNorm_yh.m` | Per-movie mesh/track alignment, cell assignment, normalized-position calculation, and diagnostics. |
| `Function/reNumCells_yh.m` | Renumbers cells consecutively across source movies. |
| `Function/getCellInfo.m` | Calculates cell area, length, estimated width, and pole boundaries. |
| `Function/getTracksInfo.m` | Collects track counts, lengths, origins, cell numbers, and flags. |
| `Function/diffAnalysis.m` | Calculates MSDs and fits diffusion/localization/power-law parameters. |
| `Function/dataToVectors.m` | Converts per-track fields into fixed-size matrices for plotting. |
| `Function/stepToVectors.m` | Calculates single-step distances, longitudinal steps, angles, normalized locations, and radius of gyration. |

## Notes

- Review the saved overlays before accepting results, especially for movies marked `badMovie` or cells not marked `filled`.
- The script assumes all tracks have enough points for the fixed summaries it creates: at least four spots for the four-frame arrays, at least twelve for the twelve-frame arrays, and at least three MSD lags for fitting.
- `autoSpotNorm.m` intentionally saves the single-day matrix output but does not create the multi-day-style `Steps` or fully analyzed `tracksFinal All` files; those are created by `combineTF.m`.
