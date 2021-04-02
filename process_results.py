import os
import numpy as np
import pandas as pd
import re


def read_and_split_s(path):
    # Reading file and storing each table as pandas dataframe
    with open(path, "r") as f:
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


def read_and_split_c(path):
    # Reading file and storing each table as pandas dataframe
    with open(path, "r") as f:
        dfs = []
        j = -1
        for line in f:
            if line[0:2] == '""':
                cols = line.strip().split(",")
                cols = ["type"] + [c.strip("\"'") for c in cols][1:]
                df = pd.DataFrame(columns=cols)
                dfs.append(df)
                j += 1
                i = 0
            elif "df" in locals():
                row = line.strip().split(",")
                row = [r.strip("\"'") for r in row]
                row = [row[0]] + [np.float(r) for r in row[1:]]
                dfs[j].loc[i] = row
                i += 1

    return dfs


def create_names(n):
    out = []
    out.append("Data")
    out.append("Distribution")
    sink = [s for s in n if "sinkhorn" in s]
    for s in sink:
        tmp = "\\(S_{"
        out.append(tmp + s[8:] + "}\\)")
    out.append("\\(MMD_u^2\\)")
    return out


def write_row(df, idx):
    val = df.loc[idx][1:].values
    row_out = str(f"{val[0]:.3f}")
    for v in val[1:]:
        row_out = row_out + " & " + f"{v:.3f}"
    return row_out + " \\\\\n"


def to_table(dfs, path, rownames, colnames):
    with open(path, "w") as f:
        f.write("\\begin{table}[ht]\n")
        f.write("\\centering\n")
        f.write("\\begin{tabular}{>{\\bfseries}r" + len(colnames) * "c" + "}\n")
        f.write("\t" + " & ".join(colnames) + " \\\\\n")
        f.write("\t\\midrule\n")
        for i, df in enumerate(dfs):
            f.write(
                "\t"
                + "\\multirow{2}*{"
                + rownames[i]
                + "} & Same & "
                + write_row(df, 0)
            )
            f.write("\t" + "& Different & " + write_row(df, 1))
            f.write("\t\\cline{2-" + str(len(colnames)) + "}\n")
        f.write("\\end{tabular}\n")
        f.write("\\caption{A Table}\n")
        f.write("\\end{table}")


def main():
    # Processing results from simulated data
    path = "./results/result_sim.out"
    dirs = ["location/sink_g.csv", "scale/sink.csv", "location/sink_l.csv"]
    for i, df in enumerate(read_and_split_s(path)):
        df.to_csv("./results/multivariate/sim/" + dirs[i], header=True, index=False)

    # Processing results from cancer data
    path = "./results/multivariate/cancer/"
    dfs = read_and_split_c(path + "result_cancer.out")
    rownames = ["Breast", "Prostate", "DLBCL"]
    colnames = [
        "Data",
        "Distribution",
        "\\(pW_2\\)",
        "\\(S_{0.1}\\)",
        "\\(S_{1}\\)",
        "\\(S_{10}\\)",
        "\\(S_{100}\\)",
        "\\(MMD_u^2\\)",
    ]
    to_table(dfs, path + "result_cancer.txt", rownames, colnames)

    # Processing univariate results from misc. data
    path = "./results/multivariate/misc/"
    dfs = read_and_split_c(path + "result_misc.out")
    rownames = ["Wines", "Ionosphere", "Iris"]
    colnames = [
        "Data",
        "Distribution",
        "\\(pW_2\\)",
        "\\(pSW_2\\)",
        "\\(S_{0.1}\\)",
        "\\(S_{1}\\)",
        "\\(S_{10}\\)",
        "\\(S_{100}\\)",
        "\\(MMD_u^2\\)",
    ]
    to_table(dfs, path + "result_misc.txt", rownames, colnames)

    # Processing multivariate results from misc. data
    path = "./results/univariate/misc/"
    dfs = read_and_split_c(path + "result_misc.out")
    rownames = ["Wines", "Ionosphere", "Iris"]
    colnames = [
        "Data",
        "Distribution",
        "\\(KS\\)",
        "\\(RT_2\\)",
        "\\(W_2\\)",
        "\\(pW_2\\)",
        "\\(S_{0.1}\\)",
        "\\(S_{1}\\)",
        "\\(S_{10}\\)",
        "\\(S_{100}\\)",
        "\\(MMD_u^2\\)",
    ]
    to_table(dfs, path + "result_misc.txt", rownames, colnames)


if __name__ == "__main__":
    main()
