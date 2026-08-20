gcd <- function(n, m) {
  stopifnot(is.numeric(n), is.numeric(m))
  if (length(n) != 1 || floor(n) != ceiling(n) ||
    length(m) != 1 || floor(m) != ceiling(m)) {
    stop("Arguments 'n', 'm' must be integer scalars.")
  }
  if (n == 0 && m == 0) {
    return(0)
  }

  n <- abs(n)
  m <- abs(m)
  if (m > n) {
    t <- n
    n <- m
    m <- t
  }
  while (m > 0) {
    t <- n
    n <- m
    m <- t %% m
  }
  return(n)
}

lcm <- function(n, m) {
  stopifnot(is.numeric(n), is.numeric(m))
  if (length(n) != 1 || floor(n) != ceiling(n) ||
    length(m) != 1 || floor(m) != ceiling(m)) {
    stop("Arguments 'n', 'm' must be integer scalars.")
  }
  if (n == 0 && m == 0) {
    return(0)
  }

  return(n / gcd(n, m) * m)
}

plcm <- function(x) {
  stopifnot(is.numeric(x))
  # if (any(floor(x) != ceiling(x)) || length(x) < 2)
  #   stop("Argument 'x' must be an integer vector of length >= 2.")

  x <- x[x != 0]
  n <- length(x)
  if (n == 0) {
    l <- 0
  } else if (n == 1) {
    l <- x
  } else if (n == 2) {
    l <- lcm(x[1], x[2])
  } else {
    l <- lcm(x[1], x[2])
    for (i in 3:n) {
      l <- lcm(l, x[i])
    }
  }
  return(l)
}
