"""
- Extract PERF lines from '.code'.
- Build a pandas dataframe
"""

import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

expr = re.compile(r"PERF +(\w+)" + r" +([0-9\.]+)" * 3, re.M)
today = time.strftime("%Y-%m-%d")


def extract_one(filename: Path) -> pd.DataFrame:
    """Extract PERF lines from a '.code' file.

    Arguments:
        filename (Path|str): Code file.

    Returns:
        DataFrame: Data with columns: (test name, command, elapsed time).
    """
    name = Path(filename).stem
    cols = ["test", "command", "cpu", "sys", "elapsed"]
    rows = []
    try:
        with open(filename) as fobj:
            content = fobj.read()
    except UnicodeDecodeError:
        print(f"ERROR: can not read {filename}, please check the testcase!")
        return pd.DataFrame()
    # icmd = 0
    cmd_counter = {}
    for mat in expr.finditer(content):
        data = list(mat.groups())
        cmd = data.pop(0)
        cmd_counter.setdefault(cmd, 0)
        cmd_counter[cmd] += 1
        # use global or per command counter
        # icmd += 1
        icmd = cmd_counter[cmd]
        data = [float(value) for value in data]
        values = name, f"{cmd}_{icmd}", *data
        rows.append(dict(zip(cols, values)))
    return pd.DataFrame(rows)


def extract_perf(resdir: Path) -> pd.DataFrame:
    """Extract PERF lines from all '.code' files from a directory.

    Arguments:
        resdir (Path|str): Results directory containing the '.code' files.

    Returns:
        DataFrame: Data with columns: (test name, command, elapsed time).
    """
    ldf = []
    for path in Path(resdir).glob("*.code"):
        ldf.append(extract_one(path))
    return pd.concat(ldf, ignore_index=True)


def compare(df1: pd.DataFrame, df2: pd.DataFrame, key: str = "cpu") -> pd.DataFrame:
    """Compare the execution times from 2 dataframes, add a 'diff' column."""
    cols = ["test", "command", key]
    dk1 = df1[cols]
    dk2 = df2[cols]
    # only keep commands, not grouped values
    dk1 = dk1[dk1.command.str.match("[A-Z0-9_]+_[0-9]+")]
    dk2 = dk2[dk2.command.str.match("[A-Z0-9_]+_[0-9]+")]
    # to restore commands order after merging
    dk1["order"] = pd.RangeIndex(len(dk1))
    dk2["order"] = pd.RangeIndex(len(dk2))

    mrg = pd.merge(dk1, dk2, on=["test", "command"], how="outer", suffixes=("_1", "_2"))
    v1 = mrg[key + "_1"]
    v2 = mrg[key + "_2"]
    mrg["diff"] = (v2 - v1) / v1 * 100.0
    mrg = mrg.sort_values(by=["order_1", "order_2"])
    cols = ["test", "command", key + "_1", key + "_2", "diff"]
    return mrg[cols]


HEADER = """
On imprime au plus {maxcmd} commandes qui se sont dégradées (ou améliorées) de plus de {perc}%
entre les dates {date1} et {date2}.
On ne compare que les commandes qui consomment plus de {tmin} secondes.
Les temps comparés sont les temps '{key}'.

Nombre de tests : {nbtest1} avant, {nbtest2} après.
Nombre de commandes comparées : {nbcmd} ({nbcmdtot} sans filtre)

"""


class PerfReport:
    """Write the report comparing to execution.

    Arguments:
        key (str): Column to be used for comparison: "cpu", "sys" or "elapsed".
    """

    _redate = re.compile("([0-9]{4}-[0-9]{2}-[0-9]{2})")

    def __init__(self, key: str = "cpu"):
        self.df1 = None
        self.df2 = None
        self.date1 = ""
        self.date2 = ""
        self.key = key
        # results
        self.cmp = None
        self.filtered = None
        self.faster = None
        self.slower = None

    def set_data(self, df1: pd.DataFrame, df2: pd.DataFrame, date1: str, date2: str):
        """Assign the dataframes to be compared.

        Arguments:
            df1 (DataFrame): First dataframe.
            df2 (DataFrame): Second dataframe.
            date1 (str): Date of the first dataframe.
            date2 (str): Date of the first dataframe.
        """
        self.df1 = df1
        self.df2 = df2
        self.date1 = date1
        self.date2 = date2
        self.cmp = compare(df1, df2, key=self.key)

    def from_csv(self, file1: Path, file2: Path):
        """Read data from csv files.

        It tries to extract the date from the filename.

        Arguments:
            file1 (Path): csv file containing the first dataframe.
            file2 (Path): csv file containing the second dataframe.
        """
        df1 = pd.read_csv(file1, index_col=0)
        df2 = pd.read_csv(file2, index_col=0)
        mat = self._redate.search(str(file1))
        if mat:
            date1 = mat.group(1)
        mat = self._redate.search(str(file2))
        if mat:
            date2 = mat.group(1)
        self.set_data(df1, df2, date1, date2)

    def show_tests(self):
        """Show the removed/added testcases."""
        # new/removed commands or testcases
        names1 = self.df1["test"].unique()
        names2 = self.df2["test"].unique()
        removed = set(names1).difference(names2)
        new = set(names2).difference(names1)
        print(f"Removed testcases: {', '.join(removed)}")
        print(f"New testcases: {', '.join(new)}")

    def show_commands(self):
        """Show the removed/added commands."""
        cmp = self.cmp
        k1, k2 = self.key + "_1", self.key + "_2"
        cmdrm = cmp[cmp[k2].isna()]
        cmdrm = cmdrm[["test", "command", k1]].reset_index(drop=True)
        cmdadd = cmp[cmp[k1].isna()]
        cmdadd = cmdadd[["test", "command", k2]].reset_index(drop=True)
        print("Removed commands:")
        print(cmdrm)
        print("Added commands:")
        print(cmdadd)

    def show_difference(self, maxcmd: int = 30, tmin: float = 1.0, perc: float = 10):
        """Show the slower and faster commands."""
        cmp = self.cmp
        k1, k2 = self.key + "_1", self.key + "_2"
        # to remove NaN and inf
        nodiff = cmp["diff"].replace([np.inf, -np.inf], np.nan).notna()
        cmp = cmp[nodiff]

        tmin = 1.0
        filt = (cmp[k1] > tmin) & (cmp[k2] > tmin)
        cmp1 = cmp[filt]
        self.filtered = cmp1

        self.faster = cmp1[cmp1["diff"] < 0.0].sort_values(by="diff", ascending=True)
        self.slower = cmp1[cmp1["diff"] > 0.0].sort_values(by="diff", ascending=False)

        infos = {}
        infos["maxcmd"] = maxcmd
        infos["perc"] = perc
        infos["tmin"] = tmin
        infos["key"] = self.key
        infos["nbtest1"] = len(self.df1["test"].unique())
        infos["nbtest2"] = len(self.df2["test"].unique())
        infos["nbcmd"] = len(self.filtered)
        infos["nbcmdtot"] = len(self.df2)
        infos["date1"] = self.date1
        infos["date2"] = self.date2
        print(HEADER.format(**infos))

        fmt_t = "{{diff:8.0f}}% {{test:<12s}} {{command:<25s}} {{{k1}:8.2f}} {{{k2}:8.2f}}"
        fmt = fmt_t.format(k1=k1, k2=k2)
        cols = ["diff", "test", "command", k1, k2]
        print("Dégradations :")
        df = self.slower[:maxcmd][cols]
        for line in df.itertuples():
            data = {}
            for col in cols:
                data[col] = getattr(line, col)
            print(fmt.format(**data))

        print("Améliorations :")
        df = self.faster[:maxcmd][cols]
        for line in df.itertuples():
            data = {}
            for col in cols:
                data[col] = getattr(line, col)
            print(fmt.format(**data))


if __name__ == "__main__":
    pd.set_option("display.max_rows", 50)

    # src = "debian-12"
    # df = extract_perf(f"/local00/tmp/perf/{src}")
    # df.to_csv(f"2026-01-14_{src}.csv")
    # df = extract_perf("/local00/tmp/perf/deb2")
    # df.to_csv("2026-01-16_debian-12.csv")

    report = PerfReport()
    report.from_csv("2026-01-14_debian-12.csv", "2026-01-16_debian-12.csv")
    report.show_difference()
