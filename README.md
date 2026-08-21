# Geodesic difference-in-differences

This repository contains the implementation and reproducibility materials for
**“Geodesic difference-in-differences”** by Yidong Zhou, Daisuke Kurisu,
Taisuke Otsu, and Hans-Georg Müller.

The repository reproduces Figures 2–7 and the numerical results reported for
the two real-data applications. The simulation results used to generate Figure
2 are included, so reproducing the paper figures does not require rerunning the
time-consuming simulations.

## Quick start

Install the required R packages:

```r
install.packages(c(
  "doSNOW", "fdadensity", "foreach", "frechet", "ggplot2", "ggtern",
  "ks", "manifold", "MASS", "patchwork", "png", "pracma", "purrr",
  "readxl", "truncnorm"
))
```

Then run all six figure scripts:

```sh
Rscript scripts/run_all_figures.R
```

The PDFs are written to `output/figures/`. Each figure can also be generated
independently:

| Paper figure | Script | Input |
|---|---|---|
| Figure 2 | `figures/figure2_simulations.R` | `em.RData`, `emmv.RData`, `en.RData` |
| Figure 3 | `figures/figure3_mortality_distributions.R` | `mortality_results.RData` |
| Figure 4 | `figures/figure4_mortality_gatt.R` | `mortality_results.RData` |
| Figure 5 | `figures/figure5_mortality_pretrend.R` | `mortality_results.RData` |
| Figure 6 | `figures/figure6_energy_sphere.R` | `figure6_treated.png`, `figure6_control.png` |
| Figure 7 | `figures/figure7_energy_ternary.R` | `energy_results_1995_2020.RData` |

The Figure 6 script combines the two three-dimensional sphere plots used in
the paper, preserving their original viewing perspective; the underlying
composition estimates are produced by `analysis/energy.R`.

## Repository structure

- `R/`: geodesic DID, geodesic transport, and utility functions.
- `simulations/`: Monte Carlo generation scripts. Full runs use 500 repetitions
  and can take substantial time.
- `analysis/`: real-data analysis scripts.
- `figures/`: one self-contained script for each paper figure, Figures 2–7.
- `data/raw/`: the public-domain EIA workbook used in the paper.
- `data/derived/`: simulation results, processed analysis results, and Figure 6
  panels.
- `output/figures/`: reproduced figures.

## Real-data analyses

### Electricity generation

The EIA workbook used in the analysis is included. Run the main analysis and
the pre-treatment comparison with:

```sh
Rscript analysis/energy.R
GDID_START_YEAR=1990 GDID_END_YEAR=1995 Rscript analysis/energy.R
```

### Age-at-death distributions

The repository includes the processed analysis results used by Figures 3–5. To
recreate the result file from the Human Mortality Database (HMD), download the
female and male period life tables at one-year ages and five-year periods,
accept the HMD user agreement, and set the two source folders:

```sh
HMD_FEMALE_DIR=/path/to/fltper_1x5 \
HMD_MALE_DIR=/path/to/mltper_1x5 \
Rscript analysis/mortality.R
```

The script preserves the paper’s treatment/control groups, smoothing settings,
and shorter available pre-periods for Greece and Slovenia. It combines the
official East and West German files in memory. Raw HMD life-table files are not
mirrored here; see [`data/README.md`](data/README.md) for provenance and
licensing.

## Simulations

The three `500 x 3` error matrices in `data/derived/` contain the simulation
results produced by `simulations/mea.R`, `simulations/mvmea.R`, and
`simulations/net.R`, respectively, and are used to generate Figure 2. Running
the simulation scripts writes new results to `output/simulations/` without
overwriting these included files. The univariate generator represents each
distribution on a fixed 101-point quantile grid.
For a short example run, use:

```sh
GDID_REPETITIONS=2 GDID_WORKERS=2 Rscript simulations/net.R
```

## License

Code is released under the MIT License. Data and derived objects retain their
source-specific terms; see [`data/README.md`](data/README.md).
