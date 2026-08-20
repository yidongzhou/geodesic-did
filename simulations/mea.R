##########################
# Gaussian distributions #
##########################
library(truncnorm)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(root, "R", "gdd.R"))
source(file.path(root, "R", "lcm.R"))
library(doSNOW)
workers <- as.integer(Sys.getenv("GDID_WORKERS", unset = "2"))
cl <- makeCluster(workers)
registerDoSNOW(cl)
clusterEvalQ(cl, library(truncnorm))
parallel::clusterSetRNGStream(cl, iseed = 1)
progress <- function(q) {
  if(q %% 10 == 0) {
    cat(sprintf("%d runs are complete\n", q))
  }
}

Q <- as.integer(Sys.getenv("GDID_REPETITIONS", unset = "500"))
nVec <- c(50, 200, 1000)
pd <- 0.25# probability of being treated
alpha1 <- alpha2 <- alpha3 <- beta <- 1
M <- 100# number of quantile intervals; each distribution uses M + 1 grid points

nu00 <- c(-3, qtruncnorm(1:(M - 1) / M, a = -3, b = 3, sd = alpha1), 3)
nu01 <- c(-3, qtruncnorm(1:(M - 1) / M, a = -3, b = 3, mean = alpha2, sd = alpha1), 3)
nu10 <- c(-3, qtruncnorm(1:(M - 1) / M, a = -3, b = 3, sd = alpha1), 3)
nu11 <- c(-3, qtruncnorm(1:(M - 1) / M, a = -3, b = 3, mean = alpha2, sd = alpha1 + beta), 3)
nu11p <- gtm(nu00, nu01, nu10, optns = list(type = 'measure', lower = -3, upper = 3))

set.seed(1)
em <- foreach(n = nVec, .combine = cbind) %:%
  foreach(icount(Q), .combine = c, .options.snow = list(progress = progress)) %dopar% {
    y00 <- y01 <- y10 <- y11 <- list()
    for(i in 1:n){
      Di <- rbinom(1, 1, pd)
      if(Di) {# treated
        y10[[length(y10) + 1]] <- c(-3, qtruncnorm(1:(M - 1) / M, a = -3, b = 3,
                                             mean = rnorm(1, mean = 0, sd = 1),
                                             sd = alpha1), 3)
        y11[[length(y11) + 1]] <- c(-3, qtruncnorm(1:(M - 1) / M, a = -3, b = 3,
                                             mean = rnorm(1, mean = alpha2, sd = 1),
                                             sd = alpha1 + beta), 3)
      } else {# control
        y00[[length(y00) + 1]] <- c(-3, qtruncnorm(1:(M - 1) / M, a = -3, b = 3,
                                             mean = rnorm(1, mean = 0, sd = 1),
                                             sd = alpha1), 3)
        y01[[length(y01) + 1]] <- c(-3, qtruncnorm(1:(M - 1) / M, a = -3, b = 3,
                                             mean = rnorm(1, mean = alpha2, sd = 1),
                                             sd = alpha1), 3)
      }
    }

    gatt <- gdd(y00, y01, y10, y11, optns = list(type = 'measure', lower = -3, upper = 3))
    sqrt(sum((nu11 - gtm(gatt$gdd$start, gatt$gdd$end, nu11p, optns = list(type = 'measure', lower = -3, upper = 3)))^2) / (M + 1))
    # sqrt(sum((nu00 - brct(y00, optns = list(type = 'measure')))^2) / (M + 1))# 0.0716 0.0510 0.0449
  }
stopCluster(cl)
dir.create(file.path(root, "output", "simulations"), recursive = TRUE, showWarnings = FALSE)
save(em, file = file.path(root, "output", "simulations", "em.RData"))

aem <- colMeans(em)
# apply(em, 2, sd)
fitm <- lm(log(aem)~log(nVec))
print(summary(fitm))# archived result used in the paper: -0.412

# ggplot(data = data.frame(x = rep(0:M / M, 4),
#                          y = c(nu00, nu01, nu10, nu11),
#                          group = rep(c("nu00", "nu01", "nu10", "nu11"), each = M + 1)),
#        aes(x = x, y = y, color = group)) +
#   geom_line(linewidth = 1) +  # Plot lines with different colors
#   labs(
#     title = "Quantile Functions",
#     x = "Quantiles",
#     y = "Values",
#     color = "Function"  # Legend title
#   ) +
#   theme_minimal() +       # Apply a clean theme
#   scale_color_manual(
#     values = c("blue", "red", "green", "purple")  # Specify colors for the curves
#   )
