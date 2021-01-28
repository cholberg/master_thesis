# SCRIPT FOR GENERATING THE SIMULATED DATA USED FOR TESTING
# ---------------------------------------------------------
library(MASS)
SEED = 24
set.seed(SEED)
GLOBAL_PATH = "/Users/christianholberg/Documents/ETH Documents/Thesis/simulations/main/data/"


# Multivariate Location - fixed n 
# P = Gaussian centered at (0, 0,...), Id variance
# Q = Gaussian centered at (log(d), 0,...), Id variance 
path <- paste0(GLOBAL_PATH, "multivariate/perm/location/dim")
NSIM <- 100
n <- 300
dseq <- seq(2, 70, by=2)
for (d in dseq) {
  filename <- paste0(path, sprintf("%02d", d), "/")
  dir.create(filename, showWarnings = FALSE)
  for (k in 1:NSIM) {
    X <- mvrnorm(n, numeric(d), diag(d))
    Y <- mvrnorm(n, c(log(d), numeric(d-1)), diag(d))
    dat <- data.frame(rbind(X, Y))
    write.csv(dat, paste0(filename, "num", k, ".csv"), row.names=FALSE)
  }
}


# Multivariate Scale - fixed n
# P = Gaussian centered at (0, 0,...), Id variance
# Q = Gaussian centered at (0, 0,...), diag((4, 1, 1,...)) variance
path <- paste0(GLOBAL_PATH, "multivariate/deviation/scale/dim")
NSIM <- 100
n <- 300
dseq <- seq(2, 42, by=2)
for (d in dseq) {
  filename <- paste0(path, sprintf("%02d", d), "/")
  dir.create(filename, showWarnings = FALSE)
  for (k in 1:NSIM) {
    X <- mvrnorm(n, numeric(d), diag(d))
    Y <- mvrnorm(n, numeric(d), diag(c(10*log(d), rep(1, d-1))))
    dat <- data.frame(rbind(X, Y))
    write.csv(dat, paste0(filename, "num", k, ".csv"), row.names=FALSE)
  }
}