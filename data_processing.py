import os
import pandas as pd
import numpy as np

GLOBAL_PATH = "/Users/christianholberg/Documents/ETH Documents/Thesis/simulations/main/data/cancer/"


def read_data(filename):
    path_read = GLOBAL_PATH + filename + ".txt"
    dat = dict()
    with open(path_read) as f:
        lines = f.readlines()
        for line in lines[:-1]:
            l = str.split(line)
            dat[l[0]] = list(map(np.float, l[1:]))
        line = str.split(lines[-1])[1:]
        a = np.unique(line)[0]
        dat["label"] = [int(l == a) for l in line]
    return pd.DataFrame(dat)


def write_data(filename):
    path_write = GLOBAL_PATH + filename + ".csv"
    dat = read_data(filename)
    dat.to_csv(path_write, header=True, index=False)


def main():
    names = ["breast", "prostate", "dlbcl"]
    for n in names:
        write_data(n)


if __name__ == "__main__":
    main()


np.unique([1, 2, 2])