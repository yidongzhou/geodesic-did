#########################################
# Networks with latent block structure. #
#########################################
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(root, "R", "gdd.R"))
library(doSNOW)
workers <- as.integer(Sys.getenv("GDID_WORKERS", unset = "2"))
cl <- makeCluster(workers)
registerDoSNOW(cl)
parallel::clusterSetRNGStream(cl, iseed = 1)
progress <- function(q) {
  if(q %% 10 == 0){
    cat(sprintf("%d runs are complete\n", q))
  }
}

Q <- as.integer(Sys.getenv("GDID_REPETITIONS", unset = "500"))
nVec <- c(50, 200, 1000)
pd <- 0.25# probability of being treated
m <- c(5, 5)
theta <- c(0.5, 0.2, 0.5)
l <- c(m[1] * (m[1] - 1) / 2, m[1] * m[2], m[2] * (m[2] - 1) / 2)
alpha1 <- alpha2 <- alpha3 <- beta <- 1

gl <- function(bk1Vec, bk12Vec, bk2Vec) {
  temp1 <- matrix(0, nrow = m[1], ncol = m[1])
  temp1[lower.tri(temp1)] <- bk1Vec
  temp1 <- temp1 + t(temp1)
  temp12 <- matrix(bk12Vec, nrow = m[1])
  temp2 <- matrix(0, nrow = m[2], ncol = m[2])
  temp2[lower.tri(temp2)] <- bk2Vec
  temp2 <- temp2 + t(temp2)
  temp <- rbind(cbind(temp1, temp12), cbind(t(temp12), temp2))
  diag(temp) <- -colSums(temp)
  temp
}

bk1Vec <- rep(-theta[1] * alpha1, l[1])
bk12Vec <- rep(-theta[2] * alpha1, l[2])
bk2Vec <- rep(-theta[3] * alpha1, l[3])
nu00 <- gl(bk1Vec, bk12Vec, bk2Vec)
bk1Vec <- rep(-theta[1] * (alpha1 + alpha2), l[1])
bk12Vec <- rep(-theta[2] * (alpha1 + alpha2), l[2])
bk2Vec <- rep(-theta[3] * (alpha1 + alpha2), l[3])
nu01 <- gl(bk1Vec, bk12Vec, bk2Vec)
bk1Vec <- rep(-theta[1] * (alpha1 + alpha3), l[1])
bk12Vec <- rep(-theta[2] * (alpha1 + alpha3), l[2])
bk2Vec <- rep(-theta[3] * (alpha1 + alpha3), l[3])
nu10 <- gl(bk1Vec, bk12Vec, bk2Vec)
bk1Vec <- rep(-theta[1] * (alpha1 + alpha2 + alpha3 + beta), l[1])
bk12Vec <- rep(-theta[2] * (alpha1 + alpha2 + alpha3 + beta), l[2])
bk2Vec <- rep(-theta[3] * (alpha1 + alpha2 + alpha3 + beta), l[3])
nu11 <- gl(bk1Vec, bk12Vec, bk2Vec)
nu11p <- gtm(nu00, nu01, nu10, optns = list(type = 'network'))

set.seed(1)
en <- foreach(n = nVec, .combine = cbind) %:%
  foreach(icount(Q), .combine = c, .options.snow = list(progress = progress)) %dopar% {
    L00 <- L01 <- L10 <- L11 <- list()
    for(i in 1:n){
      Di <- rbinom(1, 1, pd)
      if(Di) {# treated
        bk1Vec <- -rbinom(l[1], 1, theta[1]) * (alpha1 + alpha3 + runif(l[1], min = -1, max = 1))
        bk12Vec <- -rbinom(l[2], 1, theta[2]) * (alpha1 + alpha3 + runif(l[2], min = -1, max = 1))
        bk2Vec <- -rbinom(l[3], 1, theta[3]) * (alpha1 + alpha3 + runif(l[3], min = -1, max = 1))
        L10[[length(L10) + 1]] <- gl(bk1Vec, bk12Vec, bk2Vec)

        bk1Vec <- -rbinom(l[1], 1, theta[1]) * (alpha1 + alpha2 + alpha3 + beta + runif(l[1], min = -1, max = 1))
        bk12Vec <- -rbinom(l[2], 1, theta[2]) * (alpha1 + alpha2 + alpha3 + beta + runif(l[2], min = -1, max = 1))
        bk2Vec <- -rbinom(l[3], 1, theta[3]) * (alpha1 + alpha2 + alpha3 + beta + runif(l[3], min = -1, max = 1))
        L11[[length(L11) + 1]] <- gl(bk1Vec, bk12Vec, bk2Vec)
      } else {# control
        bk1Vec <- -rbinom(l[1], 1, theta[1]) * (alpha1 + runif(l[1], min = -1, max = 1))
        bk12Vec <- -rbinom(l[2], 1, theta[2]) * (alpha1 + runif(l[2], min = -1, max = 1))
        bk2Vec <- -rbinom(l[3], 1, theta[3]) * (alpha1 + runif(l[3], min = -1, max = 1))
        L00[[length(L00) + 1]] <- gl(bk1Vec, bk12Vec, bk2Vec)

        bk1Vec <- -rbinom(l[1], 1, theta[1]) * (alpha1 + alpha2 + runif(l[1], min = -1, max = 1))
        bk12Vec <- -rbinom(l[2], 1, theta[2]) * (alpha1 + alpha2 + runif(l[2], min = -1, max = 1))
        bk2Vec <- -rbinom(l[3], 1, theta[3]) * (alpha1 + alpha2 + runif(l[3], min = -1, max = 1))
        L01[[length(L01) + 1]] <- gl(bk1Vec, bk12Vec, bk2Vec)
      }
    }

    gatt <- gdd(L00, L01, L10, L11, optns = list(type = 'network'))
    sqrt(sum((nu11 - gtm(gatt$gdd$start, gatt$gdd$end, nu11p, optns = list(type = 'network')))^2))
  }
stopCluster(cl)
dir.create(file.path(root, "output", "simulations"), recursive = TRUE, showWarnings = FALSE)
save(en, file = file.path(root, "output", "simulations", "en.RData"))

aen <- colMeans(en)
fitn <- lm(log(aen)~log(nVec))
cat(sprintf("Estimated log-log slope for Figure 2: %.6f\n", coef(fitn)[2]))
