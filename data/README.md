# Data provenance and licensing

## Simulation results

`derived/em.RData`, `derived/emmv.RData`, and `derived/en.RData` contain the
simulation results produced by `simulations/mea.R`, `simulations/mvmea.R`, and
`simulations/net.R`, respectively, and are used to generate Figure 2. Each file
contains a `500 x 3` matrix of estimation errors for sample sizes 50, 200, and
1000. The files are included so that Figure 2 can be reproduced without
rerunning the full simulations.

## Human Mortality Database

`derived/mortality_results.RData` contains processed analysis results created
from the Human Mortality Database (HMD) period life tables at one-year ages and
five-year periods. It contains only the smoothed curves and geodesic DID results
needed to reproduce Figures 3–5; it does not contain the original HMD life-table
files.

HMD-produced data are available under the Creative Commons Attribution 4.0
International License. HMD asks users to download current data directly rather
than redistribute copies, because its estimates are updated. See the
[HMD user agreement](https://www.mortality.org/Data/UserAgreement) and
[citation guidelines](https://www.mortality.org/Research/CitationGuidelines).

Source attribution: Human Mortality Database, Max Planck Institute for
Demographic Research (Germany), University of California, Berkeley (U.S.A.),
and French Institute for Demographic Studies (France), available at
<https://www.mortality.org/>. The paper's data were downloaded on September 13,
2024.

The paper uses 25 countries. `analysis/mortality.R` constructs Germany by
adding the `dx` columns of the official `DEUTE` and `DEUTW` life tables within
period and age, matching the analysis input used for the paper.

## U.S. Energy Information Administration

`raw/annual_generation_state2023.xls` is the state-level electricity generation
workbook used in the analysis. It covers 1990–2023 and was released in
October 2024. The workbook is from the
[U.S. Energy Information Administration](https://www.eia.gov/electricity/data/state/).
EIA data produced by the U.S. government are in the public domain; EIA requests
source acknowledgment. See its [reuse policy](https://www.eia.gov/about/copyrights_reuse.php).

The two derived energy objects reproduce the 1995–2020 analysis and the
1990–1995 pre-treatment comparison.

`derived/figure6_treated.png` and `derived/figure6_control.png` are the saved
three-dimensional renderings used for Figure 6. The figure script combines
these panels without changing their viewing perspective.
