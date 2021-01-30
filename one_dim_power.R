# SIMULATION EXPERIMENTS: NON PARAMETRIC TWO-SAMPLE TESTING IN THE 1-DIMENSIONAL CASE
# -------------------------------------------------------------------------------------
library(parallel)
library(boot)
set.seed(24)
GLOBAL_PATH <- "/Users/christianholberg/Documents/ETH Documents/Thesis/simulations/main"

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



# TEST-STATISTICS
# Functions that calculate the test-statiscs we want to compare.
# They take as input two samples X and Y (assumed to be of equal size for simplicity)
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
  (n/2) * (sum(abs(xs-ys)^2)^(1/p) / n)
}




# ESTIMATING POWER
# We check the power of the proposed tests for simulated data consisting of normal
# samples of different mean and variance. Significance level is fixed at .05

# Example 1 - Normal, different combinations of variances and means
NUM_SIM <- 200
alpha <- .05
qsup <- qrbb(1-alpha, p="sup")
qtwo <- qrbb(1-alpha, p="two")
nseq <- c(10, 30, 100, 300)
museq <- c(0, 1, 3, 4)
sigseq <- c(1, 2, 4, 8)
theta <- list()
for (mu in museq){
  for (sig in sigseq) {
    theta[[paste0("mu", mu, "sig", sig)]] <- c(mu, sig)
  }
}
# Function for parellelization
perform_test <- function(par) {
  mu <- par[1]
  sig <- par[2]
  res_tmp <- as.data.frame(matrix(numeric(length(nseq)*3), ncol=3))
  colnames(res_tmp) <- c("Kolmogorov-Smirnov", "Ramdas-Trillos", "2-Wasserstein")
  rownames(res_tmp) <- nseq
  for (n in nseq) {
    # Drawing samples and transforming
    X <- matrix(rnorm(n*NUM_SIM), nrow=NUM_SIM)
    Y <- matrix(rnorm(n*NUM_SIM, mu, sig), nrow=NUM_SIM)
    X <- as.data.frame(t(X))
    Y <- as.data.frame(t(Y))
    # Calculating statistics
    KS <- mapply(kolmogorov_smir, X, Y)
    RT <- mapply(ramdas_trill, X, Y, MoreArgs=list(p="two"))
    W <- mapply(wasserstein, X, Y)
    # Calculating accpetance regions
    acc_KS <- qsup
    acc_RT <- qtwo
    acc_W <- mapply(qboot, X, Y, MoreArgs=list(q=1-alpha, stat=wasserstein))
    # Performing tests
    rn <- paste0(n)
    res_tmp[rn, "Kolmogorov-Smirnov"] <- sum(KS >= acc_KS) / NUM_SIM
    res_tmp[rn, "Ramdas-Trillos"] <- sum(RT >= acc_RT) / NUM_SIM
    res_tmp[rn, "2-Wasserstein"] <- sum(W >= acc_W) /NUM_SIM
  }
  res_tmp
}
# Running tests
res <- mclapply(theta, perform_test, mc.cores=8)
# Writing results to .csv files
path <- paste0(GLOBAL_PATH, "/results/univariate/sim/")
for (name in names(res)) {
  write.csv(res[[name]], paste0(path, name, ".csv"))
}