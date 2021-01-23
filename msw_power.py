import sys
import numpy as np
import pandas as pd
import torch
import gsw
import os
from multiprocessing import Pool


def unitt(dat, map="tan"):
    xnorm = (torch.sum(dat ** 2, dim=1)) ** (1 / 2)
    if map == "tan":
        a = 2 * torch.atan(xnorm) / (np.pi * xnorm)
        return torch.transpose(a * torch.transpose(dat, 0, 1), 0, 1)


def acc_region(x, n, alpha=0.05):
    R = torch.mean(torch.sum(x ** 2, dim=1)).item()
    return np.sqrt(np.log(2 / alpha) / n) + 2 * np.sqrt(2 * R / n)


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
    msw = list(
        map(
            lambda x: stat.max_gsw(x[:n, :], x[-n:, :], num_start=100, lr=0.001),
            dat,
        )
    )
    c = list(map(lambda x: acc_region(unitt(x), n), dat))
    return np.sum(np.greater(msw, c)) / NUMsim


def read_data(path):
    _, files, _ = next(os.walk(path))
    files = sorted(files)
    dirs = [path + d + "/" for d in files]
    return dirs, files


def main():
    NUMcores = 10

    # Fixed n
    # P = Gaussian centered at (0, 0,...), Id variance
    # Q = Gaussian centered at (log(d), 0,...), Id variance
    print("Running test: LOCATION, FIXED SIZE")
    path = "./data/multivariate/deviation/location/"
    dirs, files = read_data(path)
    with Pool(NUMcores) as p:
        res_msw = p.map(test_wr, dirs)
    res_msw = pd.DataFrame({"res_msw": res_msw}, index=files)
    res_msw.to_csv("./results/multivariate/sim/location/msw.csv")
    print("Done")

    # Fixed n
    # P = Gaussian centered at (0, 0,...), Id variance
    # Q = Gaussian centered at (0, 0,...), diag(10*log(d), 1,...) variance
    print("Running test: SCALE, FIXED SIZE")
    path = "./data/multivariate/deviation/scale/"
    dirs, files = read_data(path)
    with Pool(NUMcores) as p:
        res_msw = p.map(test_wr, dirs)
    res_msw = pd.DataFrame({"res_msw": res_msw}, index=files)
    res_msw.to_csv("./results/multivariate/sim/scale/msw.csv")
    print("Done")


if __name__ == "__main__":
    main()
