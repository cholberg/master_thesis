# NOTE THAT THIS CODE IS INTEDED TO RUN ON THE THE ETH EULER CLUSTER FOR 
# SCIENTIFIC COMPUTING. AS SUCH, IT IS NOT SUITED FOR LOCAL EXECUTION. 
# SINCE PERMUTATION TESTS INVOLVE A HIGH NUMBER OF REPETITIVE COMPUTATIONS, 
# BELOW CODE IS RIPE WITH OPPURTUNITIES FOR PARALLELIZATION.
# --------------------------------------------------------------------------
library(parallel)
library(LaplacesDemon)
library(kernlab)
library(boot)
library(lpSolve)
SEED <- 24
set.seed(SEED)



# Funciton for calculating cost matrix for any two empirical measures
# Input:    x - sample 1, a vector of size n
#           y - sample 2, a vector of size m
#           c - cost function
# Output: n x m cost matrix
cost_matrix <- function(x, y, c) {
  cvec <- Vectorize(function(i, j) c(x[j, ], y[i, ]))
  outer(seq_len(nrow(y)), seq_len(nrow(x)), FUN=cvec)
}


# Function for calculating entropic OT cost between any empirical measures 
# using the Sinkhorn algorithm.
# Input:    C - n x m cost matrix
#           eps - regularization parameter, real number greater than 0
#           del - convergence criteria, real number greater than 0
# Output:   optimal cost, a real number greater than 0
ent_ot <- function(C, eps=1, del=0.01) {
  n <- nrow(C)
  m <- ncol(C)
  f <- numeric(n)
  g <- numeric(m)
  ot <- sum(f)/n + sum(g)/m
  ot_prev <- 1
  d <- abs(ot-ot_prev)
  while (d > del) {
    ot_prev <- ot
    f <- sapply(seq(n), function(i) {
      -eps*log(sum(exp(log(1/m) + g/eps - C[i,]/eps)))
    })
    g <- sapply(seq(m), function(j) {
      -eps*log(sum(exp(log(1/n) + f/eps - C[,j]/eps)/n))
    })
    ot <- sum(f)/n + sum(g)/m
    d <- abs(ot-ot_prev)
  }
  ot
}


# Function for calculating entropic OT cost between any empirical measure and 
# itself using the Sinkhorn algorithm.
# Input:    C - n x n cost matrix
#           eps - regularization parameter, real number greater than 0
#           del - convergence criteria, real number greater than 0
# Output:   optimal cost, a real number greater than 0
ent_ot_id <- function(C, eps=1, del=0.01) {
  n <- nrow(C)
  f <- numeric(n)
  ot <- 2 * sum(f)/n
  ot_prev <- 1
  d <- abs(ot - ot_prev)
  while (d > del) {
    ot_prev <- ot
    f <- sapply(seq(n), function(i) {
      (f[i] - eps*log(sum(exp(f/eps - C[i,]/eps)/n)))/2
    })
    ot <- 2 * sum(f)/n
    d <- abs(ot - ot_prev)
  }
  ot
}


# Function for calculating Sinkhorn divergence between any two empirical measures
# itself using the Sinkhorn algorithm.
# Input:    C - n x m cost matrix
#           Cx - n x n cost matrix
#           Cy - m x m cost matrix
#           eps - regularization parameter, real number greater than 0
#           del - convergence criteria, real number greater than 0
# Output:   sinkhorn divergence, a real number greater than 0
sinkhorn <- function(C, Cx, Cy, eps=1, del=0.01) {
  ent_ot(C, eps, del) - 
    (ent_ot_id(Cx, eps, del) + ent_ot_id(Cy, eps, del))/2
}

# Function for calculating Wasserstein dist. between any two empirical measures
# Input:    C - n x m cost matrix
#           p - order of distance
# Output:   wasserstein distance, a real number greater than 0
wasser <- function(C) {
  # solving opt. assignment problem
  lp_sol <- lp.assign(C)
  lp_obj <- lp_sol$objval
  # computing optimal cost
  (1/nrow(C)) * lp_obj
}


# One dimensional projected 2-Wasserstein distance
#   Input:  x - sample, n x d matrix
#           y - sample, n x d matrix
#           theta - projection, vector in R^d
#   Output: a scalar greater than 0
pw <- function(x, y, theta) {
  # projecting samples + sorting
  x_s <- sort(x %*% theta)
  y_s <- sort(y %*% theta)
  # computing Wasserstein distance
  sum(abs(x_s-y_s)^2) / length(x_s)
}


# Function for sampling from (d-1)-dim unit sphere
#   Input:  n - number of samples
#           d - dimension
#   Output: n x d matrix where each row is on the unit sphere
runis <- function(n, d) {
  out <- matrix(numeric(n*d), nrow=n)
  for (i in 1:n) {
    obs <- rnorm(d)
    out[i,] <- obs / sqrt(sum(obs^2))
  }
  out
}

# Sliced 2-Wasserstein distance
sliced_wasser <- function(x, y, theta) {
  L <- NROW(theta)
  sum(apply(theta, 1, pw, x=x, y=y)) / L
}

# Funtion for performing permutation test using sinkhorn divergence
# Input:    x - sample 1, vector of size n
#           y - sample 2, vector of size m
#           alpha - significance level, number between 0 and 1
#           eps - regularization parameter, real number greater than 0
#           del - convergence criteria, real number greater than 0
#           N - number of permutations
# Output:   S3 class with the following attributes
#             res - result of test, either 0 or 1
#             stat - value of test statistic
#             quant - quantile approximated by permutations
sinkperm <- function(C, alph=.05, eps=1, del=0.01, N=1000) {
  n <- NROW(C) / 2
  m <- n
  s <- sinkhorn(
    C[1:n, (n+1):(n+m)], 
    C[1:n, 1:n], 
    C[(n+1):(n+m), (n+1):(n+m)], 
    eps, 
    del
  )
  perm <- lapply(numeric(N), function(i) sample(n+m, replace=FALSE))
  sstar <- unlist(lapply(perm, function(p) {
    sinkhorn(
      C[p[1:n], p[(n+1):(n+m)]], 
      C[p[1:n], p[1:n]], 
      C[p[(n+1):(n+m)], p[(n+1):(n+m)]],
      eps,
      del
    )
  }))
  q <- quantile(c(sstar, s), 1-alph)
  r <- as.numeric(s >= q)
  out <- list(res=r, stat=s, quant=q)
  class(out) <- "SinkhornPermutationTest"
  out
}


# Funtion for performing permutation test using sinkhorn divergence in one dimension
# Input:    x - sample 1, vector of size n
#           y - sample 2, vector of size m
#           alpha - significance level, number between 0 and 1
#           eps - regularization parameter, real number greater than 0
#           del - convergence criteria, real number greater than 0
#           N - number of permutations
# Output:   S3 class with the following attributes
#             res - result of test, either 0 or 1
#             stat - value of test statistic
#             quant - quantile approximated by permutations
sinkperm_uni <- function(x, y, alph=.05, eps=1, del=0.01, N=1000) {
  n <- NROW(x)
  m <- NROW(y)
  C <- outer(c(x, y), c(x, y), function(t, s) abs(t - s)^2)
  C <- C / max(C) # for stabilizing purposes (does not affect result)
  s <- sinkhorn(
    C[1:n, (n+1):(n+m)], 
    C[1:n, 1:n], 
    C[(n+1):(n+m), (n+1):(n+m)], 
    eps, 
    del
  )
  perm <- lapply(numeric(N), function(i) sample(n+m, replace=FALSE))
  sstar <- unlist(lapply(perm, function(p) {
    sinkhorn(
      C[p[1:n], p[(n+1):(n+m)]], 
      C[p[1:n], p[1:n]], 
      C[p[(n+1):(n+m)], p[(n+1):(n+m)]],
      eps,
      del
    )
  }))
  q <- quantile(c(sstar, s), 1-alph)
  r <- as.numeric(s >= q)
  out <- list(res=r, stat=s, quant=q)
  class(out) <- "SinkhornPermutationTest"
  out
}

# Funtion for performing permutation test using 2-Wasserstein dist.
# Input:    C - cost matrix of size n x n
#           alpha - significance level, number between 0 and 1
#           N - number of permutations
# Output:   S3 class with the following attributes
#             res - result of test, either 0 or 1
#             stat - value of test statistic
#             quant - quantile approximated by permutations
wassperm <- function(C, alph=.05, N=1000) {
  n <- NROW(C) / 2
  m <- n
  s <- wasser(C[1:n, (n+1):(n+m)])
  perm <- lapply(numeric(N), function(i) sample(n+m, replace=FALSE))
  sstar <- unlist(lapply(perm, function(p) {
    wasser(C[p[1:n], p[(n+1):(n+m)]])
  }))
  q <- quantile(c(sstar, s), 1-alph)
  r <- as.numeric(s >= q)
  out <- list(res=r, stat=s, quant=q)
  class(out) <- "WassersteinPermutationTest"
  out
}

# Funtion for performing permutation test using 2-Wasserstein dist. in one dimension
# Input:    C - cost matrix of size n x n
#           alpha - significance level, number between 0 and 1
#           N - number of permutations
# Output:   S3 class with the following attributes
#             res - result of test, either 0 or 1
#             stat - value of test statistic
#             quant - quantile approximated by permutations
wassperm_uni <- function(x, y, alph=.05, N=1000) {
  C <- outer(c(x, y), c(x, y), function(t, s) abs(t - s)^2)
  C <- C / max(C) # for stabilizing purposes (does not affect result)
  n <- NROW(C) / 2
  m <- n
  s <- wasser(C[1:n, (n+1):(n+m)])
  perm <- lapply(numeric(N), function(i) sample(n+m, replace=FALSE))
  sstar <- unlist(lapply(perm, function(p) {
    wasser(C[p[1:n], p[(n+1):(n+m)]])
  }))
  q <- quantile(c(sstar, s), 1-alph)
  r <- as.numeric(s >= q)
  out <- list(res=r, stat=s, quant=q)
  class(out) <- "WassersteinPermutationTest"
  out
}

# Funtion for performing permutation test using sliced 2-Wasserstein dist.
# Input:    x - n x d matrix
#           y - n x d matrix
#           L - number of projections
#           alpha - significance level, number between 0 and 1
#           N - number of permutations
# Output:   S3 class with the following attributes
#             res - result of test, either 0 or 1
#             stat - value of test statistic
#             quant - quantile approximated by permutations
swperm <- function(x, y, L=1000, alph=.05, N=1000) {
  n <- NROW(x)
  m <- NROW(y)
  d <- NCOL(x)
  theta <- runis(L, d)
  s <- sliced_wasser(x, y, theta)
  z <- rbind(x, y)
  perm <- lapply(numeric(N), function(i) sample(n + m, replace=FALSE))
  sstar <- unlist(lapply(perm, function(p) {
    sliced_wasser(z[p[1:n],], z[p[(n+1):(n+m)], ], theta)
  }))
  q <- quantile(c(sstar, s), 1-alph)
  r <- as.numeric(s >= q)
  out <- list(res=r, stat=s, quant=q)
  class(out) <- "SlicedWassersteinPermutationTest"
  out
}


# Funtion for performing permutation test using MMD
# Input:    x - sample 1, vector of size n
#           y - sample 2, vector of size m
#           alpha - significance level, number between 0 and 1
#           N - number of permutations
# Output:   S3 class with the following attributes
#             res - result of test, either 0 or 1
#             stat - value of test statistic
#             quant - quantile approximated by permutations
mmdperm <- function(x, y, ker="rbfdot", alph=.05, N=1000) {
  n <- NROW(x)
  m <- NROW(y)
  xm <- as.matrix(x)
  ym <- as.matrix(y)
  invisible(capture.output(
    t <- mmdstats(kmmd(xm, ym, kernel=ker))[2]
    ))
  perm <- lapply(numeric(N), function(x) sample(n+m, replace=FALSE))
  z <- rbind(xm, ym)
  tstar <- unlist(lapply(perm, function(p) {
    mmdstats(kmmd(z[p[1:n],], z[p[(n+1):(n+m)],], kernel=ker))[2]
  }))
  q <- quantile(c(tstar, t), 1-alph)
  r <- as.numeric(t >= q)
  out <- list(res=r, stat=t, quant=q)
  class(out) <- "MMDPermutationTest"
  out
}


# Funtion for performing permutation test using MMD in one dimension
# Input:    x - sample 1, vector of size n
#           y - sample 2, vector of size m
#           alpha - significance level, number between 0 and 1
#           N - number of permutations
# Output:   S3 class with the following attributes
#             res - result of test, either 0 or 1
#             stat - value of test statistic
#             quant - quantile approximated by permutations
mmdperm_uni <- function(x, y, alph=.05, N=1000) {
  n <- NROW(x)
  m <- NROW(y)
  invisible(capture.output(t <- mmdstats(kmmd(matrix(x), matrix(y)))[2]))
  perm <- lapply(numeric(N), function(x) sample(n+m, replace=FALSE))
  z <- c(x, y)
  tstar <- unlist(lapply(perm, function(p) {
    mmdstats(kmmd(matrix(z[p[1:n]]), matrix(z[p[(n+1):(n+m)]])))[2]
  }))
  tstar <- tstar[!is.na(tstar)]
  q <- quantile(c(tstar, t), 1-alph)
  r <- as.numeric(t >= q)
  out <- list(res=r, stat=t, quant=q)
  class(out) <- "MMDPermutationTest"
  out
}


# Sampling from sup_[0,1] B_t - where B is a brownain bridge on [0,1]
#   Input:  n - number of samples to draw
#           h - number of points at which B is evaluated
#   Output: a vector of length n
rbbsup <- function(n, h=1000) {
  t <- seq(0, 1, length.out = h)
  dw <- matrix(rnorm(h*n), nrow=n)/sqrt(h)
  w <- apply(dw, 1, cumsum)
  b <- abs(w - outer(t, w[h, ], "*"))
  apply(b, 2, max)
}

# Sampling from int_0^1 (B_t)^2 - where B is a brownain bridge on [0,1]
#   Input:  n - number of samples to draw
#           h - number of points at which B is evaluated
#   Output: a vector of length n
rbbell <- function(n, h=1000) {
  t <- seq(0, 1, length.out = h)
  dw <- matrix(rnorm(h*n), nrow=n)/sqrt(h)
  w <- apply(dw, 1, cumsum)
  b <- abs(w - outer(t, w[h, ], "*"))
  apply(b^2, 2, sum) / h
}

# Calculating empirical quantile for ||B_t||^p_L^p - where B is a brownian 
# bridge on [0,1] and p is either 2 or infinity
#   Input:  q - probability, number between 0 and 1
#           p - "two" for L^2 norm, "sup" for L^infty norm
#           nsim - number of simulations
#   Output: a scalar greater than 0
qrbb <- function(q, p="two", nsim=1e+5) {
  if (p=="two") {
    B <- rbbell(nsim)
  } else if (p=="sup") {
    B <- rbbsup(nsim)
  }
  quantile(B, q)
}

# Calculating quantiles for the asymptotic distribution of the test statistic
# based on the 2-Wasserstein distance using the bootstrap.
#   Input:  q - probability, number between 0 and 1
#           x - sample 1, vector of size n
#           y - sample 2, vector of size m
#           stat - test statistic, function of two samples
#           nsim - number of bootstrap simulations
#   Output: a scalar greater than 0
qboot <- function(q, x, y, stat, nsim=1000) {
  n <- length(x)
  m <- length(y)
  stat_wr <- function(dat, indices) {
    stat(dat[indices[1:n]], dat[indices[(n+1):(n+m)]])
  }
  quantile(boot(c(x, y), stat_wr, nsim)$t, q)
}


# Kolmogorov-Smirnov
#   Input:  x - sample, vector in R^n
#           y - sample, vector in R^n
#   Output: a scalar greater than 0
kolmogorov_smir <- function(x, y) {
  n <- length(x)
  z <- c(x, y)
  Fn <- Vectorize(function(t) sum(x <= t)/n)
  Gn <- Vectorize(function(t) sum(y <= t)/n)
  sqrt(n/2)*max(abs(Fn(z) - Gn(z)))
}

# Ramdas-Trillos
#   Input:  x - sample, vector in R^n
#           y - sample, vector in R^n
#           p - "two" for L^2 norm, "sup" for L^infty norm
#   Output: a scalar greater than 0
ramdas_trill <- function(x, y, p="two") {
  n <- length(x)
  xs <- sort(x)
  ys <- sort(y)
  Fn <- Vectorize(function(t) sum(x <= t)/n)
  Gn <- Vectorize(function(t) sum(y <= t)/n)
  Fq <- Vectorize(function(p) xs[which(Fn(xs) >= p)[1]])
  h <- max(10*n, 1000)
  t <- seq(0, 1, length.out=h)
  integrand <- Gn(Fq(t)) - t
  if (p == "two") {(n/2)*sum(integrand^2) / h} else {sqrt(n/2)*max(abs(integrand))}
}

# p-Wasserstein distance
#   Input:  x - sample, vector in R^n
#           y - sample, vector in R^n
#           p - order of the distance, number in [1, infty)
#   Output: a scalar greater than 0
wasserstein1d <- function(x, y, p=2) {
  n <- length(x)
  xs <- sort(x)
  ys <- sort(y)
  (n/2) * (sum(abs(xs-ys)^2)^(1/p) / n)
}

# Function for processing the results in the third simulation experiment below.
# Input:    res - a list of vectors of the size k
# Output:   a 2 x (k-1) dataframe
process_results_d <- function(res) {
  k <- length(res[[1]])
  res_mat <- matrix(unlist(res), ncol=k, byrow=TRUE)
  res_same <- res_mat[res_mat[,k]==1, 1:(k-1)]
  res_same <- res_same / NROW(res_same)
  res_diff <- res_mat[res_mat[,k]==0, 1:(k-1)]
  res_diff <- res_diff / NROW(res_diff)
  res_out <- rbind(apply(res_same, 2, sum), apply(res_diff, 2, sum))
  colnames(res_out) <- c(
    "wasserstein",
    "sinkhorn0.1", 
    "sinkhorn1", 
    "sinkhorn10", 
    "sinkhorn100",
    "mmd"
  )
  rownames(res_out) <- c("same", "different")
  res_out
}

# Function for processing the results in the second experiment below.
# Input:    res - a list of vectors of the size k
# Output:   a 2 x (k-1) dataframe
process_results_uni <- function(res) {
  res_ <- res[!(lapply(res, is.null) == TRUE)]
  N <- length(res_)
  res_ <- res_[sample(N, 200)]
  k <- length(res_[[1]])
  res_mat <- matrix(unlist(res_), ncol=k, byrow=TRUE)
  res_same <- res_mat[res_mat[,k]==1, 1:(k-1)]
  res_same <- res_same / NROW(res_same)
  res_diff <- res_mat[res_mat[,k]==0, 1:(k-1)]
  res_diff <- res_diff / NROW(res_diff)
  res_out <- rbind(apply(res_same, 2, sum), apply(res_diff, 2, sum))
  colnames(res_out) <- c(
    "kolmogorov-smirnov",
    "ramdas-trillos",
    "wasserstein",
    "pwasserstein",
    "sinkhorn0.1", 
    "sinkhorn1", 
    "sinkhorn10", 
    "sinkhorn100",
    "mmd"
  )
  rownames(res_out) <- c("same", "different")
  res_out
}


# Function for processing the results in the second experiment below.
# Input:    res - a list of vectors of the size k
# Output:   a 2 x (k-1) dataframe
process_results_mul <- function(res) {
  res_ <- res[!(lapply(res, is.null) == TRUE)]
  N <- length(res_)
  res_ <- res_[sample(N, 200)]
  k <- length(res_[[1]])
  res_mat <- matrix(unlist(res_), ncol=k, byrow=TRUE)
  res_same <- res_mat[res_mat[,k]==1, 1:(k-1)]
  res_same <- res_same / NROW(res_same)
  res_diff <- res_mat[res_mat[,k]==0, 1:(k-1)]
  res_diff <- res_diff / NROW(res_diff)
  res_out <- rbind(apply(res_same, 2, sum), apply(res_diff, 2, sum))
  colnames(res_out) <- c(
    "wasserstein",
    "sliced_wasserstein",
    "sinkhorn0.1", 
    "sinkhorn1", 
    "sinkhorn10", 
    "sinkhorn100",
    "mmd"
  )
  rownames(res_out) <- c("same", "different")
  res_out
}

# DATA INTEGRATION
# We test the usefulness in data integration tasks for 3 different 
# tumor microarray datasets
test_wr <- function(dat) {
  a <- sample(c(0, 1), 2, replace=TRUE)
  x <- dat[dat["label"]==a[1],]
  x <- as.matrix(x[sample(nrow(x), 10),])
  y <- dat[dat["label"]==a[2],]
  y <- as.matrix(y[sample(nrow(y), 10),])
  C <- cost_matrix(
    rbind(x, y), 
    rbind(x, y),
    function(s, t) sum((s-t)^2)
  )
  C <- C / max(C) # for stabilizing purposes (does not affect result)
  c(
    wassperm(C, alph=.05, N=500)$res,
    sinkperm(C, alph=.05, 0.1, N=500)$res,
    sinkperm(C, alph=.05, 1, N=500)$res,
    sinkperm(C, alph=.05, 10, N=500)$res,
    sinkperm(C, alph=.05, 100, N=500)$res,
    mmdperm(x, y, alph=.05, N=500)$res,
    as.numeric(a[1]==a[2])
  )
}

dat_names <- c("breast", "prostate", "dlbcl")
for (n in dat_names) {
  dat <- read.csv(paste0("./", n, ".csv"))
  invisible(capture.output(
    res <- process_results_d(
      mclapply(1:300, function(i) test_wr(dat), mc.cores=40)
    )
  ))
  write.csv(res, row.names=TRUE)
}


# Another test of power this time including univariate tests
alpha <- .05
qsup <- qrbb(1-alpha, p="sup")
qtwo <- qrbb(1-alpha, p="two")

test_uni_wr <- function(dat) {
  dat[] <- sapply(dat, as.numeric) # making sure everything is numeric
  l <- NROW(unique(dat["label"]))
  a <- sample(l, 2, replace = TRUE)
  k <- sample(NCOL(dat), 1)
  datA <- as.vector(dat[dat["label"] == a[1], k])
  datB <- as.vector(dat[dat["label"] == a[2], k])
  x <- datA[sample(length(datA), 10)]
  y <- datB[sample(length(datB), 10)]
  run_test <- function() {
    res <- numeric(10)
    # Univariate tests
    KS <- kolmogorov_smir(x, y)
    RT <- ramdas_trill(x, y, p="two")
    WS <- wasserstein1d(x, y)
    acc_KS <- qsup
    acc_RT <- qtwo
    acc_WS <- qboot(1-alpha, x, y, wasserstein1d)
    res[1] <- KS >= acc_KS
    res[2] <- RT >= acc_RT
    res[3] <- WS >= acc_WS
    # Multivariate tests
    res[4] <- wassperm_uni(x, y, alph=alpha, N=500)$res
    res[5] <- sinkperm_uni(x, y, alph=alpha, 0.1, N=500)$res
    res[6] <- sinkperm_uni(x, y, alph=alpha, 1, N=500)$res
    res[7] <- sinkperm_uni(x, y, alph=alpha, 10, N=500)$res
    res[8] <- sinkperm_uni(x, y, alph=alpha, 100, N=500)$res
    res[9] <- mmdperm_uni(x, y, alph=.05, N=500)$res
    # Same attribute
    res[10] <- as.numeric(a[1] == a[2])
    res
  }
  tryCatch(run_test(), error = function(e) NULL)
}


test_mul_wr <- function(dat) {
  dat[] <- sapply(dat, as.numeric) # making sure everything is numeric
  l <- NROW(unique(dat["label"]))
  a <- sample(l, 2, replace = TRUE)
  datA <- as.matrix(dat[dat["label"] == a[1], ])
  datB <- as.matrix(dat[dat["label"] == a[2], ])
  x <- datA[sample(NROW(datA), 10), ]
  y <- datB[sample(NROW(datB), 10), ]
  run_test <- function() {
    res <- numeric(8)
    C <- cost_matrix(
      rbind(x, y), 
      rbind(x, y),
      function(x, y) sum((x-y)^2)
    )
    C <- C / max(C) # for stabilizing purposes (does not affect result)
    # Multivariate tests
    res[1] <- wassperm(C, alph=alpha, N=500)$res
    res[2] <- swperm(x, y, L=1000, alph=alpha, N=500)$res
    res[3] <- sinkperm(C, alph=alpha, 0.1, N=500)$res
    res[4] <- sinkperm(C, alph=alpha, 1, N=500)$res
    res[5] <- sinkperm(C, alph=alpha, 10, N=500)$res
    res[6] <- sinkperm(C, alph=alpha, 100, N=500)$res
    res[7] <- mmdperm(x, y, alph=.05, N=500)$res
    # Same attribute
    res[8] <- as.numeric(a[1] == a[2])
    # Output result
    res
  }
  tryCatch(run_test(), error = function(e) NULL)
}

dat_names <- c("wine_processed", "ionosphere_processed", "iris_processed")
for (n in dat_names) {
  dat <- read.csv(paste0("./", n, ".csv"))
  invisible(capture.output(
    res <- process_results_mul(
      mclapply(1:250, function(i) test_mul_wr(dat), mc.cores=40)
      )
    ))
  write.csv(res, row.names=TRUE)
}
for (n in dat_names) {
  dat <- read.csv(paste0("./", n, ".csv"))
  invisible(capture.output(
    res <- process_results_uni(
      mclapply(1:250, function(i) test_uni_wr(dat), mc.cores=40)
    )
  ))
  write.csv(res, row.names=TRUE)
}