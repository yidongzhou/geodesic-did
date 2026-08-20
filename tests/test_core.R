# Fast smoke tests for the core geodesic DID implementation.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(root, "R", "lcm.R"))
source(file.path(root, "R", "gdd.R"))

euclidean <- gdd(
  y00 = list(1, 3), y01 = list(2, 4),
  y10 = list(10, 12), y11 = list(20, 22)
)
stopifnot(isTRUE(all.equal(euclidean$gdd$start, 12)),
          isTRUE(all.equal(euclidean$gdd$end, 21)))

measure <- gdd(
  y00 = list(c(0, 1, 2, 3), c(0, 1, 2, 3)),
  y01 = list(c(1, 2, 3, 4), c(1, 2, 3, 4)),
  y10 = list(c(2, 3, 4, 5), c(2, 3, 4, 5)),
  y11 = list(c(4, 5, 6, 7), c(4, 5, 6, 7)),
  optns = list(type = "measure")
)
stopifnot(all(is.finite(measure$gdd$start)),
          all(is.finite(measure$gdd$end)))

network <- gdd(
  y00 = list(diag(2), 2 * diag(2)),
  y01 = list(2 * diag(2), 3 * diag(2)),
  y10 = list(4 * diag(2), 6 * diag(2)),
  y11 = list(8 * diag(2), 10 * diag(2)),
  optns = list(type = "network")
)
stopifnot(isTRUE(all.equal(network$gdd$start, 6 * diag(2))),
          isTRUE(all.equal(network$gdd$end, 9 * diag(2))))

composition <- gdd(
  y00 = list(c(0.8, 0.2), c(0.7, 0.3)),
  y01 = list(c(0.7, 0.3), c(0.6, 0.4)),
  y10 = list(c(0.5, 0.5), c(0.4, 0.6)),
  y11 = list(c(0.3, 0.7), c(0.2, 0.8)),
  optns = list(type = "composition")
)
stopifnot(abs(sum(composition$gdd$start) - 1) < 1e-10,
          abs(sum(composition$gdd$end) - 1) < 1e-10)

set.seed(2026)
make_sample <- function(shift) {
  matrix(rnorm(60, mean = shift, sd = 0.5), ncol = 2)
}
mvmeasure <- gdd(
  y00 = list(make_sample(0), make_sample(0.1)),
  y01 = list(make_sample(0.2), make_sample(0.3)),
  y10 = list(make_sample(0.4), make_sample(0.5)),
  y11 = list(make_sample(0.6), make_sample(0.7)),
  optns = list(type = "mvmeasure", lower = c(-4, -4), upper = c(4, 4))
)
stopifnot(abs(sum(mvmeasure$gdd$start) - 1) < 1e-10,
          abs(sum(mvmeasure$gdd$end) - 1) < 1e-10,
          abs(sum(mvmeasure$omega) - 1) < 1e-10)

cat("Core smoke tests passed.\n")
