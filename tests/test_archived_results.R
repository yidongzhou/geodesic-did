# Validate archived Monte Carlo output and real-data numerical anchors.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
load_one <- function(path, name) {
  environment <- new.env(parent = emptyenv())
  load(path, envir = environment)
  environment[[name]]
}
n_values <- c(50, 200, 1000)
em <- load_one(file.path(root, "data", "derived", "em.RData"), "em")
emmv <- load_one(file.path(root, "data", "derived", "emmv.RData"), "emmv")
en <- load_one(file.path(root, "data", "derived", "en.RData"), "en")
slopes <- vapply(list(em, emmv, en), function(x) {
  unname(coef(lm(log(colMeans(x)) ~ log(n_values)))[2])
}, numeric(1))
stopifnot(isTRUE(all.equal(slopes, c(-0.4122140, -0.4730584, -0.5091448),
                                  tolerance = 1e-6)))

mortality <- load_one(
  file.path(root, "data", "derived", "mortality_results.RData"),
  "mortality_results"
)
stopifnot(
  isTRUE(all.equal(mortality$main$Female$distance, 2.114461517, tolerance = 1e-7)),
  isTRUE(all.equal(mortality$main$Male$distance, 4.258858158, tolerance = 1e-7)),
  isTRUE(all.equal(mortality$placebo$Female$distance, 0.586718557, tolerance = 1e-7)),
  isTRUE(all.equal(mortality$placebo$Male$distance, 0.932983344, tolerance = 1e-7))
)

energy_post <- load_one(
  file.path(root, "data", "derived", "energy_results_1995_2020.RData"),
  "energy_results"
)
energy_pre <- load_one(
  file.path(root, "data", "derived", "energy_results_1990_1995.RData"),
  "energy_results"
)
stopifnot(
  isTRUE(all.equal(energy_post$distance, 0.2154441, tolerance = 1e-6)),
  isTRUE(all.equal(energy_pre$distance, 0.0160604, tolerance = 1e-6))
)

cat("Archived-result checks passed.\n")
