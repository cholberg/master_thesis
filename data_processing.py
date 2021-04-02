import os
import pandas as pd
import numpy as np

GLOBAL_PATH = (
    "/Users/christianholberg/Documents/ETH Documents/Thesis/simulations/main/data/"
)


def read_data_cancer(filename):
    path_read = GLOBAL_PATH + "cancer/" + filename + ".txt"
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


def keep_continuous(dat, thr):
    cols_to_keep = []
    dat = dat.iloc[:, 1:]
    dat.columns = [*dat.columns[:-1], "label"]
    dat["label"] = np.where(dat["label"] == dat["label"].unique()[0], 2, 1)
    for col in dat.columns:
        if dat[col].nunique() >= thr:
            cols_to_keep.append(col)
    cols_to_keep.append("label")
    return dat[cols_to_keep]


def main():
    np.random.seed(24)

    # Processing cancer data
    names = ["breast", "prostate", "dlbcl"]
    for n in names:
        dat = read_data_cancer(n)
        path_write = GLOBAL_PATH + "cancer/" + n + ".csv"
        dat.to_csv(path_write, header=True, index=False)

    # Processing attribute matching data

    # Wine data
    dat = pd.read_csv(GLOBAL_PATH + "misc/wine.data", sep=",", header=None)
    dat = dat.rename(columns={0: "label"})
    dat = dat.dropna()
    dat.to_csv(GLOBAL_PATH + "misc/wine_processed.csv", index=False)

    # Ionosphere data
    dat = pd.read_csv(GLOBAL_PATH + "misc/ionosphere.csv", header=None)
    dat = keep_continuous(dat, 200)
    dat = dat.dropna()
    dat.to_csv(GLOBAL_PATH + "misc/ionosphere_processed.csv", index=False)

    # Iris data
    dat = pd.read_csv(GLOBAL_PATH + "misc/iris.data", header=None)
    dat = dat.rename(columns={4: "label"})
    dat["label"] = dat["label"].map(
        {"Iris-virginica": 1, "Iris-versicolor": 2, "Iris-setosa": 3}
    )
    dat.to_csv(GLOBAL_PATH + "misc/iris_processed.csv", index=False)

    """
    # Census data
    dat = pd.read_csv(GLOBAL_PATH + "att_matching/census_income.data", header=None)
    dat = dat.select_dtypes(include=["number"]).dropna()
    dat = dat.sample(400)
    dat.to_csv(GLOBAL_PATH + "att_matching/census.csv")
    # Forest data
    dat = pd.read_csv(GLOBAL_PATH + "att_matching/covtype.data", header=None)
    dat = dat.iloc[:, list(range(10)) + [-1]]
    dat = dat.rename(columns={54: "type"})
    dat.to_csv(GLOBAL_PATH + "att_matching/forest.csv")
    # Internet usage data
    dat = pd.read_csv(
        GLOBAL_PATH + "att_matching/final_general.data",
        sep=" ",
        header=None,
        error_bad_lines=False,
    )
    dat = dat.loc[:, dat.apply(lambda x: x.nunique()) >= 5]
    dat = dat.dropna().sample(400)
    dat.to_csv(GLOBAL_PATH + "att_matching/internet.csv")
    """


if __name__ == "__main__":
    main()
