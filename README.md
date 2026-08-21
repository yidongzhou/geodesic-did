# Geodesic difference-in-differences

This repository provides an implementation of the method introduced in
**“Geodesic difference-in-differences”** by Yidong Zhou, Daisuke Kurisu,
Taisuke Otsu, and Hans-Georg Müller. The main function, `gdd()`, supports
Euclidean, distributional, compositional, functional, network, and symmetric
positive definite matrix outcomes. Code and data for reproducing the paper's
simulations and applications are also included.

## Quick start

The four inputs to `gdd()` are lists containing the observed outcomes in the
standard difference-in-differences group-time cells:

| | Pre-treatment | Post-treatment |
|---|---|---|
| Control | `y00` | `y01` |
| Treated | `y10` | `y11` |

Each list element is one unit's outcome.

The following example applies geodesic DID to three-part compositional
outcomes. Run it from the repository root.

```r
install.packages("manifold")  # required once
source("R/lcm.R")
source("R/gdd.R")

y00 <- list(
  c(fossil = .50, nuclear = .20, renewable = .30),
  c(fossil = .40, nuclear = .25, renewable = .35)
)
y01 <- list(
  c(fossil = .45, nuclear = .20, renewable = .35),
  c(fossil = .35, nuclear = .25, renewable = .40)
)
y10 <- list(
  c(fossil = .60, nuclear = .20, renewable = .20),
  c(fossil = .55, nuclear = .25, renewable = .20)
)
y11 <- list(
  c(fossil = .35, nuclear = .20, renewable = .45),
  c(fossil = .30, nuclear = .25, renewable = .45)
)

fit <- gdd(y00, y01, y10, y11, optns = list(type = "composition"))
round(rbind(
  counterfactual = fit$gdd$start,
  observed = fit$gdd$end
), 3)
```

`fit$gdd$start` is the estimated post-treatment Fréchet mean for the treated
group in the absence of treatment, obtained by transporting the control-group
change from the treated-group pre-treatment mean. `fit$gdd$end` is the observed
treated-group post-treatment Fréchet mean. Under the conditions in the paper,
the geodesic connecting these two objects estimates the geodesic average
treatment effect on the treated. The treated-group pre-treatment Fréchet mean
is available as `fit$omega`. For Euclidean outcomes, this construction reduces
to the usual DID contrast.

## Supported outcomes

Set `optns$type` according to the outcome representation:

| `type` | Each element of an input list |
|---|---|
| `"euclidean"` | A scalar or numeric vector; this is the default |
| `"composition"` | A nonnegative vector summing to one |
| `"function"` | Function values on a common grid |
| `"measure"` | A sample from a univariate distribution |
| `"mvmeasure"` | A matrix whose rows are multivariate observations |
| `"network"` | A graph Laplacian matrix |
| `"spd"` | A symmetric positive-definite matrix |

For `type = "mvmeasure"`, supply finite `lower` and `upper` bounds in `optns`.
For univariate distributions, support bounds may be supplied but are otherwise
inferred from the data.

## Repository structure

- `R/`: geodesic DID, geodesic transport, and utility functions.
- `simulations/`: Monte Carlo generation scripts. Full runs use 500 repetitions
  and can take substantial time.
- `analysis/`: real-data analysis scripts.
- `figures/`: six self-contained scripts and their corresponding paper figures.
- `data/raw/`: the public-domain EIA workbook used in the paper.
- `data/derived/`: simulation results, processed analysis results, and Figure 6
  panels.

## Reproducing the paper results

The paper's analyses, simulations, and figures use additional R packages:

```r
install.packages(c(
  "doSNOW", "fdadensity", "foreach", "frechet", "ggplot2", "ggtern",
  "ks", "manifold", "MASS", "patchwork", "png", "pracma", "purrr",
  "readxl", "truncnorm"
))
```

### Figures

Run a figure script directly from the repository root. Each script saves its
PDF in `figures/`.

| Paper figure | Script | Input |
|---|---|---|
| Figure 2 | `figures/figure2_simulations.R` | `em.RData`, `emmv.RData`, `en.RData` |
| Figure 3 | `figures/figure3_mortality_distributions.R` | `mortality_results.RData` |
| Figure 4 | `figures/figure4_mortality_gatt.R` | `mortality_results.RData` |
| Figure 5 | `figures/figure5_mortality_pretrend.R` | `mortality_results.RData` |
| Figure 6 | `figures/figure6_energy_sphere.R` | `figure6_treated.png`, `figure6_control.png` |
| Figure 7 | `figures/figure7_energy_ternary.R` | `energy_results_1995_2020.RData` |

For example:

```sh
Rscript figures/figure2_simulations.R
```

The Figure 6 script combines the two three-dimensional sphere plots used in
the paper, preserving their original viewing perspective; the underlying
composition estimates are produced by `analysis/energy.R`.

### Real-data analyses

#### Electricity generation

The EIA workbook used in the analysis is included. Run the main analysis and
the pre-treatment comparison with:

```sh
Rscript analysis/energy.R
GDID_START_YEAR=1990 GDID_END_YEAR=1995 Rscript analysis/energy.R
```

#### Age-at-death distributions

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

### Simulations

The three `500 x 3` error matrices in `data/derived/` contain the simulation
results produced by `simulations/mea.R`, `simulations/mvmea.R`, and
`simulations/net.R`, respectively, and are used to generate Figure 2. The
simulation scripts write new results to `simulations/output/` without
overwriting these included files. The univariate generator represents each
distribution on a fixed 101-point quantile grid.

For a short example run, use:

```sh
GDID_REPETITIONS=2 GDID_WORKERS=2 Rscript simulations/net.R
```

## License

Code is released under the MIT License. Data and derived objects retain their
source-specific terms; see [`data/README.md`](data/README.md).
