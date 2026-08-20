# Geodesic difference-in-differences

This repository contains the implementation and reproducibility materials for
**“Geodesic difference-in-differences”** by Yidong Zhou, Daisuke Kurisu,
Taisuke Otsu, and Hans-Georg Müller.

The repository reproduces Figures 2–7 and the numerical results reported for
the two real-data applications. The original Monte Carlo output used in Figure
2 is included, so reproducing the paper figures does not require rerunning the
time-consuming simulations.

## Quick start

The code was checked with R 4.5.2. Install the required packages:

```r
install.packages(c(
  "doSNOW", "fdadensity", "foreach", "frechet", "ggplot2", "ggtern",
  "ks", "manifold", "MASS", "patchwork", "plotly", "pracma", "purrr",
  "readxl", "scatterplot3d", "truncnorm"
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
| Figure 6 | `figures/figure6_energy_sphere.R` | `energy_results_1995_2020.RData` |
| Figure 7 | `figures/figure7_energy_ternary.R` | `energy_results_1995_2020.RData` |

Figure 6 is a deterministic static-PDF rendering of the spherical composition
plot used in the paper.

## Reproduced numerical results

The scripts validate the following values before saving their outputs.

| Result | Reproduced value |
|---|---:|
| Figure 2 slope, univariate distributions | -0.412214 |
| Figure 2 slope, multivariate distributions | -0.473058 |
| Figure 2 slope, networks | -0.509145 |
| Mortality GATT distance, female | 2.114462 |
| Mortality GATT distance, male | 4.258858 |
| Mortality pre-treatment distance, female | 0.586719 |
| Mortality pre-treatment distance, male | 0.932983 |
| Energy GATT distance, 1995–2020 | 0.215444 |
| Energy pre-treatment distance, 1990–1995 | 0.016060 |

The energy GATT starts at `(0.431707, 0.195107, 0.373186)` and ends at
`(0.231255, 0.255217, 0.513528)`, ordered as fossil, nuclear, and renewable.

## Repository structure

- `R/`: geodesic DID, geodesic transport, and utility functions.
- `simulations/`: Monte Carlo generation scripts. Full runs use 500 repetitions
  and can take substantial time.
- `analysis/`: real-data analysis scripts.
- `figures/`: one self-contained script for each paper figure, Figures 2–7.
- `data/raw/`: the public-domain EIA workbook used in the paper.
- `data/derived/`: archived Monte Carlo output and frozen analysis objects.
- `output/figures/`: reproduced figures.
- `tests/`: fast checks for the core implementation and archived results.

## Real-data analyses

### Electricity generation

The EIA workbook used in the analysis is included. Run the main analysis and
the pre-treatment comparison with:

```sh
Rscript analysis/energy.R
GDID_START_YEAR=1990 GDID_END_YEAR=1995 Rscript analysis/energy.R
```

### Age-at-death distributions

The repository includes the frozen derived object used by Figures 3–5. To
recreate it from the Human Mortality Database (HMD), download the female and
male period life tables at one-year ages and five-year periods, accept the HMD
user agreement, and set the two source folders:

```sh
HMD_FEMALE_DIR=/path/to/fltper_1x5 \
HMD_MALE_DIR=/path/to/mltper_1x5 \
Rscript analysis/mortality.R
```

The script preserves the paper’s treatment/control groups, smoothing settings,
and shorter available pre-periods for Greece and Slovenia. It combines the
official East and West German files in memory. Raw HMD life-table files are not mirrored here; see
[`data/README.md`](data/README.md) for provenance and licensing.

## Simulations

The three archived `500 x 3` error matrices in `data/derived/` are the exact
Monte Carlo output used for Figure 2. The generation scripts write new results
to `output/simulations/`, leaving the archived objects unchanged. They preserve
the supplied computational implementation; in particular, the univariate
generator represents each distribution on a fixed 101-point quantile grid.
For a short diagnostic run, use for example:

```sh
GDID_REPETITIONS=2 GDID_WORKERS=2 Rscript simulations/net.R
```

## Verification

```sh
Rscript tests/test_core.R
Rscript tests/test_archived_results.R
```

Checksums for the archived data are in `data/SHA256SUMS`. Package versions from
the verified run are recorded in `sessionInfo.txt`.

To verify the archived files from the repository root, run
`shasum -a 256 -c data/SHA256SUMS`.

`DESCRIPTION` is a dependency manifest for this reproducibility compendium;
the repository is not structured as an installable R package.

## License

Code is released under the MIT License. Data and derived objects retain their
source-specific terms; see [`data/README.md`](data/README.md).
