"""
This script provides these features:

- Extract execution times from the '.code' files from one or several directories:

    check_perf.py --extract --output=destination-file.csv dir1 [dir2 [...]]

- Compare two or more csv files (files will be alphabetically ordered):

    check_perf.py --compare file1.csv file2.csv [file3.csv [...]]

- Check the history for a testcase (files will be alphabetically ordered):

    check_perf.py --history --name=testname --output=destination-file.csv *.csv


  The CSV file can be loaded with:

    python3
    >>> import pandas as pd
    >>> df = pd.read_csv("destination-file.csv", index_col=0)
"""

import argparse
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


def extract_perf(resdir: list[Path]) -> pd.DataFrame:
    """Extract PERF lines from all '.code' files from a directory.

    Arguments:
        resdir (list[Path|str]): Results directory containing the '.code' files.

    Returns:
        DataFrame: Data with columns: (test name, command, elapsed time).
    """
    dirs = resdir if type(resdir) in (list, tuple) else [resdir]
    ldf = []
    for dir1 in dirs:
        for path in Path(dir1).glob("*.code"):
            ldf.append(extract_one(path))
    return pd.concat(ldf, ignore_index=True)


def compare(
    df1: pd.DataFrame, df2: pd.DataFrame, key: str = "cpu", keep: list[str] = None
) -> pd.DataFrame:
    """Compare the execution times from 2 dataframes, add a 'diff' column."""
    columns = ["test", "command"]
    if not keep:
        keep = []
    cols = columns + keep + [key]
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
    if keep:
        mapping = dict([(i + "_1", i) for i in keep])
        mrg.rename(columns=mapping, inplace=True)
    cols = columns + keep + [key + "_1", key + "_2", "diff"]
    return mrg[cols]


def without_nan(df, column):
    """Remove lines containing a NaN (or inf) in a column."""
    no_nan = df[column].replace([np.inf, -np.inf], np.nan).notna()
    return df[no_nan]


HEADER = """
On imprime au plus {maxcmd} commandes qui se sont dégradées (ou améliorées) de plus de {perc}%
{hosts_compared}entre les dates {date1} et {date2}.
{host_shown}
On ne compare que les commandes qui consomment plus de {tmin} secondes.
Les temps comparés sont les temps '{key}'.

Nombre de tests : {nbtest1} avant, {nbtest2} après.
Nombre de commandes comparées : {nbcmd} ({nbcmdtot} sans filtre)
"""


regdate = re.compile("([0-9]{4}-[0-9]{2}-[0-9]{2})")


class PerfReport:
    """Write the report comparing to execution.

    Arguments:
        key (str): Column to be used for comparison: "cpu", "sys" or "elapsed".
    """

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
        mat = regdate.search(str(file1))
        if mat:
            date1 = mat.group(1)
        mat = regdate.search(str(file2))
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

    def compute_changes(self, tmin: float, perc: float):
        """Compute the slower and faster commands."""
        cmp = self.cmp
        k1, k2 = self.key + "_1", self.key + "_2"
        # to remove NaN and inf
        cmp = without_nan(cmp, "diff")

        tmin = 1.0
        filt = (cmp[k1] > tmin) & (cmp[k2] > tmin)
        cmp1 = cmp[filt]
        self.filtered = cmp1

        self.faster = cmp1[cmp1["diff"] < -perc].sort_values(by="diff", ascending=True)
        self.slower = cmp1[cmp1["diff"] > perc].sort_values(by="diff", ascending=False)

    def header_infos(self, maxcmd, tmin, perc, host_shown="", hosts_compared=""):
        """Return informations to format the header."""
        infos = {}
        infos["host_shown"] = ""
        if host_shown:
            infos["host_shown"] = f"On affiche les valeurs mesurées sur {host_shown.capitalize()}."
        infos["hosts_compared"] = ""
        if hosts_compared:
            infos["hosts_compared"] = (
                "sur " + " et ".join([i.capitalize() for i in hosts_compared]) + " "
            )
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
        return infos

    def show_difference(self, maxcmd: int = 30, tmin: float = 1.0, perc: float = 10):
        """Show the slower and faster commands."""
        self.compute_changes(tmin, perc)

        infos = self.header_infos(maxcmd, tmin, perc)
        print(HEADER.format(**infos))
        print_changes(self.slower, "Dégradations", maxcmd, self.key)
        print_changes(self.faster, "Améliorations", maxcmd, self.key)


def get_comparison(file1, file2, tmin: float, perc: float):
    report = PerfReport()
    report.from_csv(file1, file2)
    report.compute_changes(tmin, perc)
    return report


def compare_hosts(filenames: dict, maxcmd: int = 30, tmin: float = 1.0, perc: float = 10):
    """Compare 2 dates for each host.

    It expects 2 files for each host.
    The host names are deduced from the filenames which must follow this
    form: <hostname>/<date>.csv.

    Arguments:
        files (dict): 2 filenames for each host.

    Returns:
        tuple(pd.DataFrame): DataFrame of slower and faster commands.
    """
    host1, host2 = list(filenames.keys())
    data1 = get_comparison(*filenames[host1], tmin, perc)
    data2 = get_comparison(*filenames[host2], tmin, perc)

    cols = ["test", "command", "cpu_1", "cpu_2", "diff_1"]
    faster = compare(data1.faster, data2.faster, key="diff", keep=["cpu_1", "cpu_2"])
    faster = without_nan(faster, "diff")[cols]
    faster.rename(columns={"diff_1": "diff"}, inplace=True)
    slower = compare(data1.slower, data2.slower, key="diff", keep=["cpu_1", "cpu_2"])
    slower = without_nan(slower, "diff")[cols]
    slower.rename(columns={"diff_1": "diff"}, inplace=True)

    infos = data1.header_infos(maxcmd, tmin, perc, host_shown=host1, hosts_compared=(host1, host2))
    print(HEADER.format(**infos))
    print_changes(slower, "Dégradations", maxcmd)
    print_changes(faster, "Améliorations", maxcmd)
    return slower, faster


def print_changes(df: pd.DataFrame, legend: str, maxcmd: int = 30, key: str = "cpu"):
    """Print a DataFrame containing columns: test, command, key_1, key_2, diff"""
    k1, k2 = key + "_1", key + "_2"
    fmt_t = "{{diff:8.0f}}% {{test:<12s}} {{command:<25s}} {{{k1}:8.2f}} {{{k2}:8.2f}}"
    fmt = fmt_t.format(k1=k1, k2=k2)
    cols = ["diff", "test", "command", k1, k2]
    print(f"{legend}:")
    df = df[:maxcmd][cols]
    for line in df.itertuples():
        data = {}
        for col in cols:
            data[col] = getattr(line, col)
        print(fmt.format(**data))


if __name__ == "__main__":
    pd.set_option("display.max_rows", 50)

    parser = argparse.ArgumentParser(
        usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--extract", action="store_true", help="Extract execution times from '.code' files"
    )
    parser.add_argument("--compare", action="store_true", help="Compare two CSV files")
    parser.add_argument(
        "--compare-host",
        dest="compare_host",
        action="store_true",
        help="Get host name from CSV files, expecting 4 CSV files (2 for each host), "
        "must follow this form: <hostname>/<date>.csv",
    )
    parser.add_argument("--history", action="store_true", help="History for a testcase")
    parser.add_argument(
        "-o", "--output", action="store", help="Filename of the CSV file for '--extract'"
    )
    parser.add_argument("--name", action="store", help="Name of the testcase")
    parser.add_argument(
        "args", metavar="FILES/DIRS", nargs="*", help="List of files or directories"
    )
    args = parser.parse_args()
    if args.extract:
        if not args.output:
            parser.error("please define the result file using '-o'/'--output' option")
        if not args.args:
            parser.error("at least one directory must be passed in argument")
        df = extract_perf(args.args)
        df.to_csv(args.output)
    elif args.compare:
        if len(args.args) != 2:
            parser.error("expecting 2 CSV files for comparison")
        files = sorted(args.args, reverse=True)
        file2 = files.pop(0)
        while files:
            file1, file2 = file2, files.pop(0)
            report = PerfReport()
            report.from_csv(file2, file1)
            report.show_difference()
    elif args.compare_host:
        if len(args.args) != 4:
            parser.error("expecting 4 CSV files for comparison")
        filenames = {}
        for csv in args.args:
            csv = Path(csv)
            mat = regdate.search(csv.name)
            assert mat, f"filename must match YYYY-mm-dd format, not {csv.name}"
            host = csv.parent.name
            filenames.setdefault(host, [])
            filenames[host].append(csv)
        for host, files in filenames.items():
            if len(files) != 2:
                raise ValueError(f"host '{host}': expecting 2 files, {len(files)} provided")
        slower, faster = compare_hosts(filenames)
    elif args.history:
        if not args.output:
            parser.error("please define the result file using '-o'/'--output' option")
        if not args.args:
            parser.error("at least one CSV file must be passed in argument")
        files = sorted(args.args)
        dft = None
        for csv in files:
            csv = Path(csv)
            print(f"reading {csv}...", end=" ")
            date = csv.stem
            mat = regdate.search(str(csv))
            if mat:
                date = mat.group(1)
            cols = ["command", "cpu"]
            df = pd.read_csv(csv, index_col=0)
            dfi = df[df["test"] == args.name][cols]
            dfi.rename(columns={"cpu": date}, inplace=True)
            if dft is None:
                dft = dfi
                continue
            dft = pd.merge(dft, dfi, how="outer", on="command")
            print(f"current size {dft.shape}")
        dft.to_csv(args.output)
