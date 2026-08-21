###############################
# Multivariate distributions  #
###############################

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(root, "R", "gdd.R"))
source(file.path(root, "R", "lcm.R"))
library(doSNOW)

workers <- as.integer(Sys.getenv("GDID_WORKERS", unset = "2"))
cl <- makeCluster(workers)
registerDoSNOW(cl)
parallel::clusterSetRNGStream(cl, iseed = 1)

progress <- function(q) {
  if (q %% 10 == 0) {
    cat(sprintf("%d runs are complete\n", q))
  }
}

Q <- as.integer(Sys.getenv("GDID_REPETITIONS", unset = "500"))
nVec <- c(50, 200, 1000)
pd <- 0.25  # probability of being treated
alpha2 <- 1
beta <- 0.2
Sigma <- 0.5 * diag(2)
N <- 500    # number of observations per (i, t)
nPopUnits <- 1000 # number of units per group for population targets

# helper to sample from the bivariate Gaussian mixture
rmix_bvn <- function(mu1, mu2, pi, n, Sigma) {
  z <- rbinom(n, size = 1, prob = pi)
  n1 <- sum(z == 1)
  n2 <- n - n1
  x1 <- if (n1 > 0) MASS::mvrnorm(n1, mu = mu1, Sigma = Sigma) else matrix(numeric(0), ncol = 2)
  x2 <- if (n2 > 0) MASS::mvrnorm(n2, mu = mu2, Sigma = Sigma) else matrix(numeric(0), ncol = 2)
  rbind(x1, x2)
}

# population-level targets for mvmeasure: nu00, nu01, nu10, nu11 and reference point nu11p
# Use the known Gaussian mixture densities evaluated on a fixed grid (no KDE / preprocessing).
set.seed(1)

# fixed 10x10 grid in R^2, matching optns$lower/upper used in gdd ([-4,4] in each dimension)
lower_vec <- c(-4, -4)
upper_vec <- c(4, 4)
gx <- seq(lower_vec[1], upper_vec[1], length.out = 10)
gy <- seq(lower_vec[2], upper_vec[2], length.out = 10)
grid <- as.matrix(expand.grid(gx, gy))  # 100 x 2

# closed-form bivariate normal density with Sigma = 0.5 I_2
phi_bvn <- function(x, mu) {
  # x: matrix M x 2, mu: length-2 vector
  diff <- sweep(x, 2, mu, FUN = "-")
  (1 / pi) * exp(-rowSums(diff^2))
}

dens00_list <- dens01_list <- dens10_list <- dens11_list <- vector("list", nPopUnits)
for (i in 1:nPopUnits) {
  # control units (D = 0)
  xi0 <- MASS::mvrnorm(1, mu = c(0, 0), Sigma = diag(2))
  mu1_00 <- xi0 + alpha2 * 0 * c(1, 1)
  mu2_00 <- xi0 - alpha2 * 0 * c(1, 1)
  pi_00  <- 0.5 + beta * 0 * 0
  f00 <- pi_00  * phi_bvn(grid, mu1_00) + (1 - pi_00)  * phi_bvn(grid, mu2_00)

  mu1_01 <- xi0 + alpha2 * 1 * c(1, 1)
  mu2_01 <- xi0 - alpha2 * 1 * c(1, 1)
  pi_01  <- 0.5 + beta * 0 * 1
  f01 <- pi_01  * phi_bvn(grid, mu1_01) + (1 - pi_01)  * phi_bvn(grid, mu2_01)

  # treated units (D = 1)
  xi1 <- MASS::mvrnorm(1, mu = c(0, 0), Sigma = diag(2))
  mu1_10 <- xi1 + alpha2 * 0 * c(1, 1)
  mu2_10 <- xi1 - alpha2 * 0 * c(1, 1)
  pi_10  <- 0.5 + beta * 1 * 0
  f10 <- pi_10  * phi_bvn(grid, mu1_10) + (1 - pi_10)  * phi_bvn(grid, mu2_10)

  mu1_11 <- xi1 + alpha2 * 1 * c(1, 1)
  mu2_11 <- xi1 - alpha2 * 1 * c(1, 1)
  pi_11  <- 0.5 + beta * 1 * 1
  f11 <- pi_11  * phi_bvn(grid, mu1_11) + (1 - pi_11)  * phi_bvn(grid, mu2_11)

  dens00_list[[i]] <- f00 / sum(f00)
  dens01_list[[i]] <- f01 / sum(f01)
  dens10_list[[i]] <- f10 / sum(f10)
  dens11_list[[i]] <- f11 / sum(f11)
}

nu00 <- brct(dens00_list, optns = list(type = "mvmeasure"))
nu01 <- brct(dens01_list, optns = list(type = "mvmeasure"))
nu10 <- brct(dens10_list, optns = list(type = "mvmeasure"))
nu11 <- brct(dens11_list, optns = list(type = "mvmeasure"))
nu11p <- gtm(nu00, nu01, nu10, optns = list(type = "mvmeasure"))

emmv <- foreach::foreach(n = nVec, .combine = cbind) %:%
  foreach::foreach(icount(Q), .combine = c, .options.snow = list(progress = progress)) %dopar% {
    y00 <- y01 <- y10 <- y11 <- list()

    for (i in 1:n) {
      Di <- rbinom(1, 1, pd)
      xi <- MASS::mvrnorm(1, mu = c(0, 0), Sigma = diag(2))

      # t = 0
      mu1_0 <- xi + alpha2 * 0 * c(1, 1)
      mu2_0 <- xi - alpha2 * 0 * c(1, 1)
      pi_0  <- 0.5 + beta * Di * 0
      Y_i0  <- rmix_bvn(mu1_0, mu2_0, pi_0, N, Sigma)

      # t = 1
      mu1_1 <- xi + alpha2 * 1 * c(1, 1)
      mu2_1 <- xi - alpha2 * 1 * c(1, 1)
      pi_1  <- 0.5 + beta * Di * 1
      Y_i1  <- rmix_bvn(mu1_1, mu2_1, pi_1, N, Sigma)

      if (Di) {  # treated
        y10[[length(y10) + 1]] <- Y_i0
        y11[[length(y11) + 1]] <- Y_i1
      } else {   # control
        y00[[length(y00) + 1]] <- Y_i0
        y01[[length(y01) + 1]] <- Y_i1
      }
    }

    gatt <- gdd(y00, y01, y10, y11,
                optns = list(type = "mvmeasure",
                             lower = lower_vec,
                             upper = upper_vec))

    # transport along the estimated geodesic to the population reference point nu11p
    beta_hat_ref <- gtm(gatt$gdd$start, gatt$gdd$end, nu11p, optns = list(type = "mvmeasure"))

    # error: Fisher--Rao distance (geodesic angle on the unit sphere)
    ab <- sum(sqrt(nu11) * sqrt(beta_hat_ref))
    ab <- max(min(ab, 1), -1)  # numerical guard
    acos(ab)
  }

stopCluster(cl)

dir.create(file.path(root, "simulations", "output"), recursive = TRUE, showWarnings = FALSE)
save(emmv, file = file.path(root, "simulations", "output", "emmv.RData"))

aemmv <- colMeans(emmv)
fitmmv <- lm(log(aemmv)~log(nVec))
cat(sprintf("Estimated log-log slope for Figure 2: %.6f\n", coef(fitmmv)[2]))
