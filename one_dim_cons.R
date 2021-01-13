# SIMULATION EXPERIMENTS: NON PARAMETRIC TWO-SAMPLE TESTING IN THE 1-DIMENSIONAL CASE
# -------------------------------------------------------------------------------------


# HELPER FUNCTIONS
# Some functions that will be useful throughout

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

# Calculating acceptance region for test 1-Wasserstein distance
#   Input:  n - size of samples, natural number greater than 2
#           alpha - significance level, number between 0 and 1
#   Ouput:  a scalar greater than 0
acceptance_reg <- function(n, alpha=.05) {
  nf<- floor(n/2)
  R <- choose(2*nf, nf)/(4^nf)
  2*R + sqrt(log(1/alpha)/n)
}

# Function mapping R^n onto the unit ball
#   Input:  x - a vector in R^n
#   Ouput:  a vector in R^n with norm less than 1
unit_trans <- function(x) {
  xnorm <- abs(x)
  a <- 2 * atan(xnorm) / (pi * xnorm)
  x * a
}




# TEST-STATISTICS
# We create functions that perform the calculate the test-statiscs we want to compare.
# They take as input two samples X and Y (assumed to be of equal size for simplicity),
# and potentially some hyperparameters

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
wasserstein <- function(x, y, p=2) {
  n <- length(x)
  xs <- sort(x)
  ys <- sort(y)
  sum(abs(xs-ys)^2)^(1/p) / n
}




# CHECKING CONSISTENCY
# We check the consitency of the proposed tests for two different simulated
# examples and one real world example. Significance level is fixed at .05

# Initializing
NUM_SIM <- 100
alpha <- .05
nseq <- c(10, 30, 100, 300)
qsup <- qrbb(1-alpha, p="sup")
qtwo <- qrbb(1-alpha, p="two")

# Example 1 - Standard Normal
res_norm <- as.data.frame(matrix(numeric(length(nseq)*4), ncol=4))
colnames(res_norm) <- c("Kolmogorov-Smirnov", "Ramdas-Trillos, p=infty", 
                   "Ramdas-Trillos, p=2", "Wasserstein")
rownames(res_norm) <- nseq
for (n in nseq) {
  # Drawing samples and transforming
  X <- matrix(rnorm(n*NUM_SIM), nrow=NUM_SIM)
  Y <- matrix(rnorm(n*NUM_SIM), nrow=NUM_SIM)
  Xprime <- as.data.frame(t(unit_trans(X)))
  Yprime <- as.data.frame(t(unit_trans(Y)))
  X <- as.data.frame(t(X))
  Y <- as.data.frame(t(Y))
  # Calculating statistics
  KS <- mapply(kolmogorov_smir, X, Y)
  RTsup <- mapply(ramdas_trill, X, Y, MoreArgs=list(p="sup"))
  RTtwo <- mapply(ramdas_trill, X, Y, MoreArgs=list(p="two"))
  WS <- mapply(wasserstein, Xprime, Yprime, MoreArgs=list(p=1))
  # Calculating accpetance regions
  acc_KS <- qsup
  acc_RTsup <- qsup
  acc_RTtwo <- qtwo
  acc_WS <- acceptance_reg(n, alpha)
  # Performing tests
  rn <- paste0(n)
  res_norm[rn, "Kolmogorov-Smirnov"] <- sum(KS >= acc_KS) / NUM_SIM
  res_norm[rn, "Ramdas-Trillos, p=infty"] <- sum(RTsup >= acc_RTsup) / NUM_SIM
  res_norm[rn, "Ramdas-Trillos, p=2"] <- sum(RTtwo >= acc_RTtwo) / NUM_SIM
  res_norm[rn, "Wasserstein"] <- sum(WS >= acc_WS) / NUM_SIM
}

# Example 2 - Uniform
res_unif <- as.data.frame(matrix(numeric(length(nseq)*4), ncol=4))
colnames(res_unif) <- c("Kolmogorov-Smirnov", "Ramdas-Trillos, p=infty", 
                   "Ramdas-Trillos, p=2", "Wasserstein")
rownames(res_unif) <- nseq
for (n in nseq) {
  # Drawing samples and transforming
  X <- matrix(runif(n*NUM_SIM), nrow=NUM_SIM)
  Y <- matrix(runif(n*NUM_SIM), nrow=NUM_SIM)
  X <- as.data.frame(t(X))
  Y <- as.data.frame(t(Y))
  # Calculating statistics
  KS <- mapply(kolmogorov_smir, X, Y)
  RTsup <- mapply(ramdas_trill, X, Y, MoreArgs=list(p="sup"))
  RTtwo <- mapply(ramdas_trill, X, Y, MoreArgs=list(p="two"))
  WS <- mapply(wasserstein, X, Y, MoreArgs=list(p=1))
  # Calculating accpetance regions
  acc_KS <- qsup
  acc_RTsup <- qsup
  acc_RTtwo <- qtwo
  acc_WS <- acceptance_reg(n, alpha)
  # Performing tests
  rn <- paste0(n)
  res_unif[rn, "Kolmogorov-Smirnov"] <- sum(KS >= acc_KS) / NUM_SIM
  res_unif[rn, "Ramdas-Trillos, p=infty"] <- sum(RTsup >= acc_RTsup) / NUM_SIM
  res_unif[rn, "Ramdas-Trillos, p=2"] <- sum(RTtwo >= acc_RTtwo) / NUM_SIM
  res_unif[rn, "Wasserstein"] <- sum(WS >= acc_WS) / NUM_SIM
}
