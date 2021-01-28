library(future.apply)
library(kernlab)
SEED <- 24
set.seed(SEED)
options(future.globals.maxSize= 1000*1024^2)


# TEST BASDE ON LARGE DEVIATION BOUND
# -------------------------------------------------------------------

# LOCATION

# Fixed n=300
# P = Gaussian centered at (0, 0,...), Id variance
# Q = Gaussian centered at (log(d), 0,...), Id variance
path <- "./data/multivariate/deviation/location/"
dirs <- sort(list.files(path)[grep("dim", list.files(path))])
NUMdat <- length(dirs)
NUMsim <- length(list.files(paste0(path, dirs[1])))
# Loading data
dat <- vector("list", NUMdat)
for (i in 1:NUMdat) {
  dat_tmp <- vector("list", NUMsim)
  path_read <- paste0(path, dirs[i], "/")
  files <- list.files(path_read)
  for (j in 1:NUMsim) {
    file_read <- paste0(path_read, files[j])
    dat_tmp[[j]] <- as.matrix(read.csv(file_read))
  }
  dat[[i]] <- dat_tmp
}
# Wrapper function
n <- dim(dat[[1]][[1]])[1] / 2
test_wr <- function(dat) {
  out <- lapply(dat, function(x) as.numeric(H0(kmmd(x[1:n,], x[-(1:n),]))))
  sum(unlist(out)) / NUMsim
}
# Running tests
invisible(capture.output(res_mmd <- unlist(future_lapply(dat, test_wr, future.seed=SEED))))
names(res_mmd) <- dirs
write.csv(data.frame(res_mmd), "./results/multivariate/sim/location/mmd.csv")


# SCALE

# Fixed n=400
# P = Gaussian centered at (0, 0,...), Id variance
# Q = Gaussian centered at (0, 0,...), diag(rep(20*log(d), d/2), rep(1, d/2)) variance
path <- "./data/multivariate/deviation/scale/"
dirs <- sort(list.files(path)[grep("dim", list.files(path))])
NUMdat <- length(dirs)
NUMsim <- length(list.files(paste0(path, dirs[1])))
# Loading data
dat <- vector("list", NUMdat)
for (i in 1:NUMdat) {
  dat_tmp <- vector("list", NUMsim)
  path_read <- paste0(path, dirs[i], "/")
  files <- list.files(path_read)
  for (j in 1:NUMsim) {
    file_read <- paste0(path_read, files[j])
    dat_tmp[[j]] <- as.matrix(read.csv(file_read))
  }
  dat[[i]] <- dat_tmp
}
# Wrapper function
n <- dim(dat[[1]][[1]])[1] / 2
test_wr <- function(dat) {
  out <- lapply(dat, function(x) as.numeric(H0(kmmd(x[1:n,], x[-(1:n),]))))
  sum(unlist(out)) / NUMsim
}
# Running tests
invisible(capture.output(res_mmd <- unlist(future_lapply(dat, test_wr, future.seed=SEED))))
names(res_mmd) <- dirs
write.csv(data.frame(res_mmd), "./results/multivariate/sim/scale/mmd.csv")



# TEST BASDE ON PERMUTATION
# -------------------------------------------------------------------


