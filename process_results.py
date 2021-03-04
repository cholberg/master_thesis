import os
import numpy as np
import pandas as pd

"./results/multivariate/sim/"


def read_and_split(path):
    fname = [f for f in os.listdir(path) if "result" in f][0]
    fname = path + fname

    # Reading file and storing each table as pandas dataframe
    with open(fname, "r") as f:
        dfs = []
        j = -1
        for line in f:
            if line[0:2] == '""':
                cols = line.strip().split(",")
                cols = ["dim"] + [c.strip("\"'") for c in cols][1:]
                df = pd.DataFrame(columns=cols)
                dfs.append(df)
                j += 1
                i = 0
            elif "df" in locals():
                row = line.strip().split(",")
                row = np.array([r.strip("\"'") for r in row]).astype(np.float)
                dfs[j].loc[i] = row
                i += 1

    return dfs


def main():
    path = "./results/multivariate/sim/"
    dirs = ["location/", "scale/"]
    for i, df in enumerate(read_and_split(path)):
        df.to_csv(path + dirs[i] + "sinkhorn.csv", header=True, index=False)


if __name__ == "__main__":
    main()
