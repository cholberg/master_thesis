# NOTE THAT THIS CODE IS INTEDED TO RUN ON THE THE ETH EULER CLUSTER FOR 
# SCIENTIFIC COMPUTING. AS SUCH, IT IS NOT SUITED FOR LOCAL EXECUTION. 
# SINCE PERMUTATION TESTS INVOLVE A HIGH NUMBER OF REPETITIVE COMPUTATIONS, 
# BELOW CODE IS RIPE WITH OPPURTUNITIES FOR PARALLELIZATION.
# --------------------------------------------------------------------------
library(parallel)
library(MASS)
library(kernlab)
library(transport)
library(lpSolve)
library(LaplacesDemon)
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


# Function for calculating squared Wasserstein dist. between any two empirical measures
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
sliced_wasser <- function(x, y, L=1000) {
  n <- dim(x)[1]
  d <- dim(x)[2]
  theta <- runis(L, d)
  sum(apply(theta, 1, pw, x=x, y=y)) / L
}


# Funtion for performing permutation test using sinkhorn divergence
# Input:    C - cost matrix of size 2n x 2n
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
  s <- sliced_wasser(x, y, L)
  z <- rbind(x, y)
  perm <- lapply(numeric(N), function(i) sample(n + m, replace=FALSE))
  sstar <- unlist(lapply(perm, function(p) {
    sliced_wasser(z[p[1:n],], z[p[(n+1):(n+m)], ], L)
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


# POWER: LOCATION
# Fixed n=100
# X ~ Gaussian with (0,0,...) mean and Id variance
# Y ~ Gaussian with (1,0,...) mean and Id variance
NSIM <- 100
n <- 100
dseq <- seq(2, 2002, by=100)


test_wr <- function(dseq, n) {
  dmax <- max(dseq)
  x <- mvrnorm(n, numeric(dmax), diag(dmax))
  y <- mvrnorm(n, c(1, numeric(dmax-1)), diag(dmax))
  run_test <- function(d) {
    C <- cost_matrix(
      rbind(x[, 1:d], y[, 1:d]), 
      rbind(x[, 1:d], y[, 1:d]),
      function(s, t) sum((s-t)^2)
    )
    C <- C / max(C) # for stabilizing purposes (does not affect result)
    c(
      wassperm(C, alph=.05, N=500)$res,
      swperm(x[, 1:d], y[, 1:d], L=1000, alph=.05, N=500)$res,
      sinkperm(C, alph=.05, eps=0.1, N=500)$res,
      sinkperm(C, alph=.05, eps=1, N=500)$res,
      sinkperm(C, alph=.05, eps=10, N=500)$res,
      sinkperm(C, alph=.05, eps=100, N=500)$res,
      mmdperm(x[, 1:d], y[, 1:d], alph=.05, N=500)$res
      )
  }
  matrix(
    unlist(mclapply(dseq, run_test, mc.cores=2)), 
    ncol=7, 
    byrow=TRUE,
    )
}


invisible(
  capture.output(
    res <- Reduce(
      '+', 
      mclapply(1:NSIM, function(i) test_wr(dseq, n), mc.cores=20)
      ) / NSIM
  )
)
colnames(res) <- c(
  "wasser",
  "sw",
  "sinkhorn0.1", 
  "sinkhorn1", 
  "sinkhorn10", 
  "sinkhorn100", 
  "mmd"
  )
rownames(res) <- dseq
write.csv(res, row.names=TRUE)



# POWER: SCALE
# Fixed n=100
# X ~ Gaussian with (0,0,...) mean and Id variance
# Y ~ Gaussian with (0,0,...) mean and diag(4,0,0,...) variance
NSIM <- 100
n <- 200
dseq <- seq(2, 102, by=5)


test_wr <- function(dseq, n) {
  dmax <- max(dseq)
  x <- mvrnorm(n, numeric(dmax), diag(dmax))
  y <- mvrnorm(n, numeric(dmax), diag(c(4, rep(1, dmax-1))))
  run_test <- function(d) {
    C <- cost_matrix(
      rbind(x[, 1:d], y[, 1:d]), 
      rbind(x[, 1:d], y[, 1:d]),
      function(s, t) sum((s-t)^2)
    )
    C <- C / max(C) # for stabilizing purposes (does not affect result)
    c(
      wassperm(C, alph=.05, N=500)$res,
      swperm(x[, 1:d], y[, 1:d], L=1000, alph=.05, N=500)$res,
      sinkperm(C, alph=.05, eps=0.1, N=500)$res,
      sinkperm(C, alph=.05, eps=1, N=500)$res,
      sinkperm(C, alph=.05, eps=10, N=500)$res,
      sinkperm(C, alph=.05, eps=100, N=500)$res,
      mmdperm(x[, 1:d], y[, 1:d], alph=.05, N=500)$res
    )
  }
  matrix(
    unlist(mclapply(dseq, run_test, mc.cores=2)), 
    ncol=7, 
    byrow=TRUE,
  )
}


invisible(
  capture.output(
    res <- Reduce(
      '+', 
      mclapply(1:NSIM, function(i) test_wr(dseq, n), mc.cores=20)
    ) / NSIM
  )
)
colnames(res) <- c(
  "wasser",
  "sw",
  "sinkhorn0.1", 
  "sinkhorn1", 
  "sinkhorn10", 
  "sinkhorn100", 
  "mmd"
)
rownames(res) <- dseq
write.csv(res, row.names=TRUE)


# POWER: LOCATION
# Fixed n=100
# X ~ Laplacian with (0,0,...) mean and Id scale
# Y ~ Laplacian with (1,0,...) mean and Id scale
NSIM <- 100
n <- 100
dseq <- seq(2, 2002, by=100)

test_wr <- function(dseq, n) {
  dmax <- max(dseq)
  x <- rmvl(n, numeric(dmax), diag(dmax))
  y <- rmvl(n, c(1, numeric(dmax-1)), diag(dmax))
  run_test <- function(d) {
    C <- cost_matrix(
      rbind(x[, 1:d], y[, 1:d]), 
      rbind(x[, 1:d], y[, 1:d]),
      function(s, t) sum((s-t)^2)
    )
    C <- C / max(C) # for stabilizing purposes (does not affect result)
    c(
      wassperm(C, alph=.05, N=500)$res,
      swperm(x[, 1:d], y[, 1:d], L=1000, alph=.05, N=500)$res,
      sinkperm(C, alph=.05, eps=0.1, N=500)$res,
      sinkperm(C, alph=.05, eps=1, N=500)$res,
      sinkperm(C, alph=.05, eps=10, N=500)$res,
      sinkperm(C, alph=.05, eps=100, N=500)$res,
      mmdperm(x[, 1:d], y[, 1:d], ker="laplacedot", alph=.05, N=500)$res
    )
  }
  matrix(
    unlist(mclapply(dseq, run_test, mc.cores=2)), 
    ncol=7, 
    byrow=TRUE,
  )
}


invisible(
  capture.output(
    res <- Reduce(
      '+', 
      mclapply(1:NSIM, function(i) test_wr(dseq, n), mc.cores=20)
    ) / NSIM
  )
)
colnames(res) <- c(
  "wasser",
  "sw",
  "sinkhorn0.1", 
  "sinkhorn1", 
  "sinkhorn10", 
  "sinkhorn100", 
  "mmd"
)
rownames(res) <- dseq
write.csv(res, row.names=TRUE)
