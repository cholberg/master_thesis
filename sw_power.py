import sys
import numpy as np
import pandas as pd
import torch
import gsw
import os
from scipy.special import binom
from multiprocessing import Pool


def unitt(dat, map="tan"):
    xnorm = (torch.sum(dat ** 2, dim=1)) ** (1 / 2)
    if map == "tan":
        a = 2 * torch.atan(xnorm) / (np.pi * xnorm)
        return torch.transpose(a * torch.transpose(dat, 0, 1), 0, 1)


def acc_region(n, alpha=0.05):
    nf = n // 2
    R = binom(2 * nf, nf) / (4 ** nf)
    return 2 * R + np.sqrt(np.log(1 / alpha) / n)


def test_wr(path):
    files = os.listdir(path)
    NUMsim = len(files)
    dat = []
    for f in files:
        file_read = path + "/" + f
        samp = torch.tensor(pd.read_csv(file_read).values, dtype=torch.float)
        dat.append(unitt(samp))
    n = dat[0].shape[0] // 2
    stat = gsw.GSW(use_cuda=False)
    sw = list(map(lambda x: stat.gsw(x[:n, :], x[-n:, :]), dat))
    c = acc_region(n)
    return np.sum(np.greater(sw, c)) / NUMsim


def read_data(path):
    _, dirs, _ = next(os.walk(path))
    dirs = sorted(dirs)
    dirs = [path + d + "/" for d in dirs]
    return dirs


def main():
    NUMcores = 10

    # Fixed n
    # P = Gaussian centered at (0, 0,...), Id variance
    # Q = Gaussian centered at (log(d), 0,...), Id variance
    print("Running test: LOCATION, FIXED SIZE")
    path = "./data/multivariate/deviation/location/"
    dirs = read_data(path)
    with Pool(NUMcores) as p:
        res_sw = p.map(test_wr, dirs)
    res_sw = pd.DataFrame({"res_sw": res_sw}, index=dirs)
    res_sw.to_csv("./results/multivariate/sim/location/sw.csv")
    print("Done")

    # Fixed n
    # P = Gaussian centered at (0, 0,...), Id variance
    # Q = Gaussian centered at (0, 0,...), diag(10*log(d), 1,...) variance
    print("Running test: SCALE, FIXED SIZE")
    path = "./data/multivariate/deviation/scale/"
    dirs = read_data(path)
    with Pool(NUMcores) as p:
        res_sw = p.map(test_wr, dirs)
    res_sw = pd.DataFrame({"res_sw": res_sw}, index=dirs)
    res_sw.to_csv("./results/multivariate/sim/scale/sw.csv")
    print("Done")


if __name__ == "__main__":
    main()