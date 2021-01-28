library(parallel)
library(MASS)
library(kernlab)
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
  f_prev <- rep(1, n)
  g_prev <- rep(1, m)
  dist <- sum(abs(f - f_prev)) + sum(abs(g - g_prev))
  while (dist > del) {
    f <- sapply(seq(n), function(i) -eps*log(sum(exp(g/eps - C[i,]/eps)/m)))
    g <- sapply(seq(m), function(j) -eps*log(sum(exp(f/eps - C[,j]/eps)/n)))
    dist <- sum(abs(f - f_prev)) + sum(abs(g - g_prev))
    f_prev <- f
    g_prev <- g
  }
  sum(f)/n + sum(g)/m
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
  f_prev <- rep(1, n)
  dist <- sum(abs(f - f_prev))
  while (dist > del) {
    f <- sapply(seq(n), function(i) {
      (f[i] - eps*log(sum(exp(f/eps - C[i,]/eps)/n)))/2
    })
    dist <- sum(abs(f - f_prev))
    f_prev <- f
  }
  2 * sum(f)/n
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
sinkperm <- function(x, y, alph=.05, eps=1, del=0.01, N=1000) {
  n <- NROW(x)
  m <- NROW(y)
  C <- cost_matrix(
    rbind(x, y), 
    rbind(x, y),
    function(x, y) sqrt(sum((x-y)^2))
    )
  s <- sinkhorn(
    C[1:n, (n+1):(n+m)], 
    C[1:n, 1:n], 
    C[(n+1):(n+m), (n+1):(n+m)], 
    eps, 
    del
    )^2 * (m*n/(m+n))
  perm <- lapply(numeric(N), function(x) sample(n+m, replace=FALSE))
  sstar <- unlist(lapply(perm, function(p) {
    sinkhorn(
      C[p[1:n], p[(n+1):(n+m)]], 
      C[p[1:n], p[1:n]], 
      C[p[(n+1):(n+m)], p[(n+1):(n+m)]],
      eps,
      del
      )^2 * (m*n/(m+n))
  }))
  q <- quantile(sstar, 1-alph)
  r <- as.numeric(s >= q)
  out <- list(res=r, stat=s, quant=q)
  class(out) <- "SinkhornPermutationTest"
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
mmdperm <- function(x, y, alph=.05, N=1000) {
  n <- NROW(x)
  m <- NROW(y)
  invisible(capture.output(t <- mmdstats(kmmd(x, y))[2]))
  perm <- lapply(numeric(N), function(x) sample(n+m, replace=FALSE))
  z <- rbind(x, y)
  tstar <- unlist(lapply(perm, function(p) {
    mmdstats(kmmd(z[p[1:n],], z[p[(n+1):(n+m)],]))[2]
    }))
  q <- quantile(tstar, 1-alph)
  r <- as.numeric(t >= q)
  out <- list(res=r, stat=t, quant=q)
  class(out) <- "MMDPermutationTest"
  out
}



# POWER: LOCATION
# Fixed n=200
# X ~ Gaussian with (0,0,...) mean and Id variance
# Y ~ Gaussian with (1,0,...) mean and Id variance
NSIM <- 100
n <- 200
dseq <- seq(2, 2002, by=100)


test_wr <- function(dseq, n) {
  dmax <- max(dseq)
  x <- mvrnorm(n, numeric(dmax), diag(dmax))
  y <- mvrnorm(n, c(1, numeric(dmax-1)), diag(dmax))
  run_test <- function(d) {
    c(
      sinkperm(x[, 1:d], y[, 1:d], alph=.05, 0.1, N=500)$res,
      sinkperm(x[, 1:d], y[, 1:d], alph=.05, 1, N=500)$res,
      sinkperm(x[, 1:d], y[, 1:d], alph=.05, 10, N=500)$res,
      sinkperm(x[, 1:d], y[, 1:d], alph=.05, 100, N=500)$res,
      mmdperm(x[, 1:d], y[, 1:d], alph=.05, N=500)$res
      )
  }
  matrix(
    unlist(mclapply(dseq, run_test, mc.cores=2)), 
    ncol=5, 
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
  "sinkhorn0.1", 
  "sinkhorn1", 
  "sinkhorn10", 
  "sinkhorn100", 
  "mmd"
  )
rownames(res) <- dseq
write.csv(res, row.names=TRUE)



# POWER: SCALE
# Fixed n=200
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
    c(
      sinkperm(x[, 1:d], y[, 1:d], alph=.05, 0.1, N=500)$res,
      sinkperm(x[, 1:d], y[, 1:d], alph=.05, 1, N=500)$res,
      sinkperm(x[, 1:d], y[, 1:d], alph=.05, 10, N=500)$res,
      sinkperm(x[, 1:d], y[, 1:d], alph=.05, 100, N=500)$res,
      mmdperm(x[, 1:d], y[, 1:d], alph=.05, N=500)$res
    )
  }
  matrix(
    unlist(mclapply(dseq, run_test, mc.cores=2)), 
    ncol=5, 
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
  "sinkhorn0.1", 
  "sinkhorn1", 
  "sinkhorn10", 
  "sinkhorn100", 
  "mmd"
)
rownames(res) <- dseq
write.csv(res, row.names=TRUE)