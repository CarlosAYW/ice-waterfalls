# Ice Waterfalls

> **Reviewing the master's thesis?** Start at **[CarlosAYW/waterfall-ice-thesis](https://github.com/CarlosAYW/waterfall-ice-thesis)**. It maps every thesis figure to its script and data and includes this repository as a submodule pinned to the thesis version (tag `thesis-submission`). That tag holds the model code exactly as used for the thesis runs, including an indexing error in the solar term of the ice core that was corrected afterwards on `main` (details in README 8 of the thesis repository).

Interactive icefall map and model output for a master thesis project.

## Repository structure

- `assets/`: homepage source (`index.html`) and static assets copied into the GitHub Pages build.
- `scripts/`: runtime build scripts used by GitHub Actions.
- `data/`: input, cache, and derived data required by the model/build.
  - `data/Glacier/`: glacier inventory shapefile used by the CAP calculation.
  - `data/inca_nordtirol/point_timeseries/`: cached INCA point time series used by UID model plots.
- `analysis/`: research, QA, and thesis-helper scripts that are not part of the normal Pages build.
- `site/`: generated GitHub Pages output. It is ignored except for `.gitkeep`.
- `ALT/`: local archive for old scripts, scratch files, and RStudio state. It is ignored by Git.

## GitHub Pages build

The Pages site is built by `.github/workflows/build_site.yml` on pushes to `main` and `main-test`.

The workflow runs:

1. `scripts/00_build_plots_all.R`
2. `scripts/02_build_list_page.R`
3. copy `assets/index.html` and other assets into `site/`
4. `scripts/01_build_map.R`

On `main`, the generated `site/` folder is deployed to GitHub Pages. On `main-test`, it is uploaded as a workflow artifact for checking.

## Local smoke test

From the repository root:

```bash
Rscript scripts/00_build_plots_all.R --uids=1
Rscript scripts/02_build_list_page.R
Rscript scripts/01_build_map.R
```

Then serve the generated site locally:

```bash
python -m http.server 8000 --directory site
```

## Add new icefalls

Fill one or more rows in `add_new/new_icefalls.csv`, then run from the
repository root:

```bash
Rscript scripts/add_new_icefalls.R
```

If `Rscript` is not in `PATH` on Windows, use:

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" scripts/add_new_icefalls.R
```

Leave `uid` empty for new icefalls; set an existing `uid` to update/recalculate
that icefall. The command updates the main table and the fixed derived
parameters such as height/aspect, station assignment, wind vulnerability,
topographic sun, route structure, and cold-air-pooling tables. Add `--dry-run`
to preview the affected UIDs without writing files.

## Rerun Oetztal with the 50 cm raster

The Oetztal 50 cm raster is expected at:

```text
data/DEM/DOM_Oetztal_50cm.tif
```

When this file exists, the DEM catalog uses it before the 5 m Tirol DGM.

For a quick test on one UID first, without rebuilding sun tables, CAP, or the
site:

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" scripts/add_new_icefalls.R --rerun-oetztal --uids=1 --skip-sun --skip-cap --skip-site --force-structure
```

To recalculate all existing Oetztal icefalls from the main table and rebuild
the model plots and map:

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" scripts/add_new_icefalls.R --rerun-oetztal --force-structure --run-models --build-map
```

## Notes for pushing

Do not commit generated caches, `site/`, local DEM source files, RStudio state, or local analysis outputs. The `.gitignore` keeps these out of Git.

The repository already contains large tracked data files, but no currently tracked file is above GitHub's 100 MB single-file limit.

## Copyright and use

Copyright (c) 2026. All rights reserved.

This repository was created as part of a master thesis. Code, models, scripts, texts, graphics, and other contents are protected by copyright.

Without prior written permission, copying, redistribution, publication, modification, derivative works, use in other projects or teaching material, and commercial use are not permitted.

Use is limited to reading, review, and assessment in the context of the master thesis.
