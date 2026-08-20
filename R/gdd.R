#' @title Geodesic Difference-in-Differences
#' @description Difference-in-differences for random objects lying in a geodesic space.
#' @param y00 A list of random objects from the control group observed in the pre-treatment time period.
#' @param y01 A list of random objects from the control group observed in the post-treatment time period.
#' @param y10 A list of random objects from the treated group observed in the pre-treatment time period.
#' @param y11 A list of random objects from the treated group observed in the post-treatment time period.
#' @param optns A list of options control parameters specified by
#'   \code{list(name = value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{type}{type of random objects. 'euclidean', 'composition', 'function',
#' 'measure', 'mvmeasure', 'network', and 'spd' are supported. The default is
#' 'euclidean'.}
#' \item{lower}{Lower support bound. For univariate measures, an omitted bound
#' is inferred from the inputs. Multivariate measures require finite bounds.}
#' \item{upper}{Upper support bound. For univariate measures, an omitted bound
#' is inferred from the inputs. Multivariate measures require finite bounds.}
#' }
#' @return A \code{gdd} object, which is a list containing the following components:
#' \describe{
#' \item{gdd}{geodesic difference-in-differences estimator.}
#' \item{omega}{the reference point.}
#' }

gdd <- function(y00, y01, y10, y11, optns = list()) {
  if (is.null(y00) || is.null(y01) || is.null(y10) || is.null(y11)) {
    stop("requires the input of y00, y01, y10, and y11")
  }
  if (!is.list(y00) || !is.list(y01) || !is.list(y10) || !is.list(y11)) {
    stop("y00, y01, y10, and y11 must be lists")
  }
  if (length(y00) != length(y01) || length(y10) != length(y11)) {
    stop("y00 and y01 must have the same length, and y10 and y11 must have the same length")
  }
  if (is.null(optns$type)) {
    optns$type <- 'euclidean'
  }
  if(is.null(optns$lower)) {
    optns$lower <- -Inf
  }
  if(is.null(optns$upper)) {
    optns$upper <- Inf
  }
  if(optns$type == 'measure') {
    if (!is.finite(optns$lower)) {
      optns$lower <- min(unlist(c(y00, y01, y10, y11)))
    }
    if (!is.finite(optns$upper)) {
      optns$upper <- max(unlist(c(y00, y01, y10, y11)))
    }
    y00 <- harmonize_measure(y00)
    y01 <- harmonize_measure(y01)
    y10 <- harmonize_measure(y10)
    y11 <- harmonize_measure(y11)
  }
  if(optns$type == 'mvmeasure') {
    y00 <- harmonize_mvmeasure(y00, optns)
    y01 <- harmonize_mvmeasure(y01, optns)
    y10 <- harmonize_mvmeasure(y10, optns)
    y11 <- harmonize_mvmeasure(y11, optns)
  }
  nu00 <- brct(y00, optns)
  nu01 <- brct(y01, optns)
  nu10 <- brct(y10, optns)
  nu11 <- brct(y11, optns)
  omega <- nu10# reference point

  mu1 <- gtm(nu00, nu01, omega, optns)
  mu2 <- nu11
  if(optns$type == 'composition') {
    names(mu1) <- names(y00[[1]])
    names(mu2) <- names(y11[[1]])
    names(omega) <- names(y10[[1]])
  }
  # For mvmeasure outcomes, keep the discrete densities normalized to sum to
  # one. The Fisher--Rao transport below operates on their unit-norm square
  # roots; multiplying by the grid size would violate that input contract.

  res <- list()
  res$gdd <- list(start = mu1, end = mu2)
  res$omega <- omega
  class(res) <- 'gdd'
  res
}

# Preprocess measure input: list of vectors (variable length) -> list of length-M vectors
harmonize_measure <- function(y) {
  n <- length(y)
  N <- sapply(y, length)
  if (length(unique(N)) == 1L) {
    return(lapply(y, sort))
  }
  M <- min(plcm(N), n * max(N), 5000)
  lapply(seq_len(n), function(i) {
    residual <- M %% N[i]
    if (residual) {
      sort(c(rep(y[[i]], each = M %/% N[i]), sample(y[[i]], residual)))
    } else {
      sort(rep(y[[i]], each = M %/% N[i]))
    }
  })
}

# Preprocess mvmeasure input: list of n_i x d matrices -> list of density vectors on a 10^d grid
harmonize_mvmeasure <- function(y, optns) {
  n <- length(y)
  d <- ncol(y[[1]])
  sample_sizes <- sapply(y, nrow)

  if (length(unique(sample_sizes)) == 1) {
    M <- sample_sizes[1]
    y_extended <- y
  } else {
    size_range <- max(sample_sizes) - min(sample_sizes)
    size_ratio <- size_range / min(sample_sizes)

    if (size_ratio < 0.2) {
      M <- min(sample_sizes)
      y_extended <- lapply(y, function(yi) yi[sample(nrow(yi), M), , drop = FALSE])
    } else {
      M <- min(plcm(sample_sizes), n * max(sample_sizes), 5000)
      y_extended <- lapply(seq_along(y), function(i) {
        yi <- y[[i]]
        n_i <- nrow(yi)
        residual <- M %% n_i
        if (residual) {
          indices <- c(rep(1:n_i, each = M %/% n_i), sample(1:n_i, residual))
          yi[indices, , drop = FALSE]
        } else {
          rep_indices <- rep(1:n_i, each = M %/% n_i)
          yi[rep_indices, , drop = FALSE]
        }
      })
    }
  }

  gridsize <- rep(10, d)
  lower <- optns$lower
  upper <- optns$upper
  if (length(lower) == 1L) lower <- rep(lower, d)
  if (length(upper) == 1L) upper <- rep(upper, d)
  lapply(seq_len(n), function(i) {
    fit <- ks::kde(x = y_extended[[i]], gridsize = gridsize, xmin = lower, xmax = upper)
    dens <- as.vector(fit$estimate)
    dens <- pmax(dens, 0)
    dens / sum(dens)
  })
}

# barycenter of a list of random objects
brct <- function(y, optns, w = rep(1 / length(y), length(y))) {# y is a list
  n <- length(y)
  if (optns$type %in% c('composition', 'mvmeasure')) {
    mfd <- structure(1, class = 'Sphere')
    yM <- matrix(unlist(y), ncol = n)
    brct <- c(manifold::frechetMean(mfd = mfd, X = sqrt(yM), weight = w, maxit = 1e04))^2
  } else if (optns$type %in% c('euclidean', 'function', 'measure', 'network', 'spd')) {
    brct <- purrr::reduce(lapply(1:n, function(i) w[i] * y[[i]]), `+`)
  }
  brct
}

# geodesic transport
gtm <- function(alpha, beta, omega, optns) {
  if (optns$type == 'measure') {
    if(alpha[1] > omega[1]) {
      omega[1] <- alpha[1]
    }
    if(omega[length(omega)] > alpha[length(alpha)]) {
      omega[length(omega)] <- alpha[length(omega)]
    }
    zeta <- pracma::spinterp(x = alpha, y = beta, xp = omega)
    zeta[1] <- optns$lower
    zeta[length(zeta)] <- optns$upper
  } else if (optns$type %in% c('composition', 'mvmeasure')) {
    alpha <- sqrt(alpha)
    beta <- sqrt(beta)
    omega <- sqrt(omega)
    # log--parallel--exp on the unit sphere (both use norm-standardized vectors)
    ab <- sum(alpha * beta)
    ao <- sum(alpha * omega)
    if (1 + ao < 1e-10) stop("The minimizing geodesic from alpha to omega is not unique (omega ≈ -alpha).")
    ab <- max(min(ab, 1), -1)
    theta <- acos(ab)
    s <- sqrt(max(1 - ab^2, 0))
    if (theta < 1e-12 || s < 1e-12) {
      zeta0 <- omega
    } else {
      u <- (beta - ab * alpha) / s
      v <- theta * u
      vo <- sum(v * omega)
      ptv <- v - (vo / (1 + ao)) * (alpha + omega)
      nv <- sqrt(sum(ptv^2))
      if (nv < 1e-12) {
        zeta0 <- omega
      } else {
        zeta0 <- cos(nv) * omega + sin(nv) * (ptv / nv)
      }
    }
    zeta <- pmax(zeta0, 0)
    if (sum(zeta^2) < 1e-15) {
      zeta <- numeric(length(zeta))
      zeta[which.max(zeta0)] <- 1
    } else {
      zeta <- zeta / sqrt(sum(zeta^2))
    }
    zeta <- zeta^2
  } else if (optns$type %in% c('euclidean', 'function', 'network', 'spd')) {
    zeta <- omega + beta - alpha
  }
  zeta
}
