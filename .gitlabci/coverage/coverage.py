#!/usr/bin/env python3
# coding=utf-8

"""
Check code coverage.
"""

import argparse
import os
import os.path as osp
import pickle
import re
import sys
import time
from shutil import copyfile
from pathlib import Path

from .features import Features

USAGE = """%(prog)s [options]

Check code coverage (keywords combinations and elementary calculations).
"""

EPILOG = """
Examples:
- for periodic coverage checking:
    %(prog)s --save

- compare to a specific analysis:
    %(prog)s --previous=2021-08-14

- export results as a text file similarly to legacy scripts:
    %(prog)s --savetxt
"""


# srcdir: data, src, validation, devtools
# installdir
# root: dest
class CoverageAnalysis:
    """Check code coverage."""

    MESS_EXT = ".mess.txt"

    def __init__(
        self,
        root: Path,
        srcdir: Path,
        installdir: Path,
        resdir: Path,
        previous: str,
        limit: int,
        verbose: bool = False,
        force: bool = False,
    ):
        self._root = root
        self._resdir = resdir
        self._force = force
        self._features = Features(installdir, srcdir, self.log)
        self._all = set()
        self._tested_by = {}  # feature: [(time, testcase name)]
        self._data = self.result_file("tested")
        self._previous = previous
        self._limit = limit
        self._verbose = verbose
        self._report = []

    def log(self, *args, **kwargs):
        """Logging function"""
        if self._verbose and (not kwargs.get("need_tty") or sys.stdout.isatty()):
            kwargs.pop("need_tty", None)
            print(*args, **kwargs)

    def report(self, text=""):
        """Add text to the report."""
        self._report.append(text)
        self.log(text)

    def check(self):
        """Check code coverage."""
        self.to_be_tested()
        self.parse()

    def result_file(self, typ: str, prefix: str = None):
        """Return path to a result file.

        Arguments:
            typ (str): file type
            prefix (str, optional): file name prefix, a timestamp is used by default.

        Returns:
            str: Path to file.
        """
        prefix = prefix or time.strftime("%Y-%m-%d")
        return osp.join(self._root, prefix + "." + typ)

    def _show_grouped(self, features: Features, tested_by: bool = None):
        grp = {}
        for feature in features:
            parent, base = _split(feature)
            grp.setdefault(base, {})
            grp[base].setdefault(True, [])
            grp[base].setdefault(False, [])
            tested = bool(tested_by and feature in tested_by)
            grp[base][tested].append(parent)

        single = []
        others = []
        for base, values in grp.items():
            if len(values[True]) + len(values[False]) == 1:
                parent = (values[True] + values[False])[0]
                single.append(_join(parent, base))
            else:
                others.append(base)
        for feature in sorted(single):
            self._show(feature, tested_by)

        for base in sorted(others):
            self.report()
            for is_tested in (True, False):
                label = "est testé dans" if is_tested else "n'est pas testé dans"
                if grp[base][is_tested]:
                    self.report(f"  * {base} {label} :")
                    for parent in sorted(grp[base][is_tested]):
                        self._show(parent, tested_by, level=3)

    def _search_long(self):
        test_long = {}
        for feature in self._tested_by:
            tests = sorted(self._tested_by[feature])[0]
            if tests[0] >= self._limit:
                test_long[feature] = [tests]
        return test_long

    def check_diff(self):
        """Check code coverage evolution since previous analysis."""
        last = self.result_file("tested", self._previous)
        if osp.exists(last):
            with open(last, "rb") as fpick:
                prev_all = pickle.load(fpick)
                prev_tested = pickle.load(fpick)
        else:
            print(f"WARNING: can not find previous analysis ({last})")
            return

        vers = self._features.get_version()
        new = self._all.difference(prev_all)
        now_tested = set(self._tested_by.keys())
        removed = prev_all.difference(self._all)

        if vers:
            self.report()
            self.report(" " * 12 + f"-- CODE_ASTER VERSION {vers} --")
        self.report()
        self.report("1. Anciennes fonctionnalités :")
        lst = set(prev_tested.keys()).difference(now_tested).difference(removed)
        self.report()
        self.report(f"- Anciennes fonctionnalités qui ne sont plus testées : {len(lst)}")
        # self.report("(nom du test le plus court qui vérifiait cette fonctionnalité)")
        for feature in sorted(lst):
            self._show(feature, prev_tested)

        self.report()
        self.report(f"2. Nouvelles fonctionnalités : {len(new)}")
        lst = new.difference(now_tested)
        self.report()
        self.report(f"- Nouvelles fonctionnalités qui ne sont pas testées : {len(lst)}")
        self._show_grouped(lst)

        self.report()
        lst = now_tested.difference(prev_tested.keys())
        self.report(f"- Fonctionnalités testées (et qui ne l'étaient pas) : {len(lst)}")
        for feature in sorted(lst):
            self._show(feature, self._tested_by)

        self.report()
        self.report(f"- Anciennes fonctionnalités supprimées du code : {len(removed)}")
        for feature in sorted(removed):
            self._show(feature)

        self.report()
        if self._limit:
            self.report("3. Statistiques")
            self.report()
            lst = self._search_long()
            self.report(
                f"- Fonctionnalités uniquement testées par des tests de plus "
                f"de {self._limit} secondes : {len(lst)}"
            )
            for feature in sorted(lst):
                self._show(feature, lst)

    def to_be_tested(self):
        """Build the list of all features."""
        self._all = self._features.get_all()

    def parse(self):
        """Loop on testcases."""
        self.log(f"INFO: searching test results from {self._resdir}")
        files = self._resdir.glob("*.code")
        # ignore export.mess/export.code...
        files = sorted([i for i in files if not i.name.startswith("export")])
        total = len(files)
        self.log(f"INFO: parsing {total} code files")
        for i, code in enumerate(files):
            i += 1
            test = code.stem
            self.log(f"parsing ({i}/{total}): {test:<60s}", end="\r", need_tty=True)
            sys.stdout.flush()
            self.check_one(code)
        self.log(" " * 80, end="\r", need_tty=True)
        self.log(f"INFO: {len(files)} code files parsed")

    def check_one(self, filename: Path):
        """Check a testcase coverage."""
        test = filename.with_suffix("")
        elaps = get_total_job(test.with_suffix(self.MESS_EXT))
        ref = elaps, test.name
        with open(filename, "rb") as fcode:
            content = fcode.read().decode(errors="replace")
        self._parse_code(ref, content)

    def _parse_code(self, ref, content):
        """Parse a 'code' file.

        Arguments:
            ref (str): identifier of the testcase.
            content (str): content of the 'code' file.
        """

        def _mark(feature):
            self._tested_by.setdefault(feature, [])
            if ref not in self._tested_by[feature]:
                self._tested_by[feature].append(ref)

        re_test = re.compile("^ *TEST *", re.M)
        re_calc = re.compile(r"&&CALCUL +(?P<option>\w+) +(?P<elt>\w+)", re.M)
        re_cmde = re.compile(
            r"^(?P<cmde>[^&]\w+) +(?P<fact>[-\w]+) +(?P<simp>\w+)" r" +(?P<value>.*)", re.M
        )
        content = re_test.sub("", content)
        calculs = []
        for mat in re_calc.finditer(content):
            phen = self._features.get_model(mat.group("elt"))
            calculs.extend([i + "/" + mat.group("option") for i in phen])
        for key in calculs:
            key = "te:" + key
            # assert key in self._all, f"{ref[1]}: {key}: inconsistent testcase"
            if key not in self._all:
                print(f"ERROR: {ref[1]}: {key}: inconsistent testcase", mat.groups())
                continue
            _mark(key)

        for mat in re_cmde.finditer(content):
            key0 = "kw:" + "/".join(mat.groups()[:3])
            values = mat.group("value")
            found = key0 in self._all
            if found:
                _mark(key0)
            elif not values:
                # hidden commands, indirectly called (not in Cata/Commands)
                # print(f"WARNING: unknown keywords: {key0}")
                pass
            # print(key0 in self._all, key0, mat.groups())
            if values:
                try:
                    values = eval(values)
                except Exception:
                    # unexpected values
                    continue
                if type(values) not in (list, tuple):
                    values = [values]
                for value in values:
                    key = f"{key0}#{value}"
                    if key in self._all:
                        _mark(key)

    def save(self):
        """Save current analysis."""
        with open(self.result_file("tested"), "wb") as fobj:
            self.log(f"INFO: writing {fobj.name}")
            pickle.dump(self._all, fobj)
            pickle.dump(self._tested_by, fobj)
        copyfile(self.result_file("tested"), self.result_file("tested", "last"))

        with open(self.result_file("results"), "w") as fout:
            self.log(f"INFO: writing {fout.name}")
            fout.write("\n".join(self._report))
        copyfile(self.result_file("results"), self.result_file("results", "last"))

    def savetxt(self):
        """Save text files."""
        with open(self.result_file("all_features"), "w") as fobj:
            self.log(f"INFO: writing {fobj.name}")
            fobj.write("\n".join(sorted(self._all)))

        not_tested = self._all.difference(set(self._tested_by.keys()))
        with open(self.result_file("not_tested"), "w") as fobj:
            self.log(f"INFO: writing {fobj.name}")
            fobj.write("\n".join(sorted(not_tested)))

    def _show(self, feature, tested_by=None, level=1):
        light = feature.replace("kw:", "").replace("te:", "")
        indent = "  " * level
        if tested_by and tested_by.get(feature):
            test = tested_by[feature]
            strtest = ", ".join([f"{name}~{elaps:.1f}s" for elaps, name in sorted(test)[:4]])
            self.report(f"{indent}{light} ({strtest})")
        else:
            self.report(f"{indent}{light}")


def get_total_job(filename):
    """Extract the execution time from a 'mess' file.

    Arguments:
        filename (str): Filename.

    Returns:
        float: The execution time (or 0.).
    """
    elaps = 0.0
    if osp.exists(filename):
        expr = r"TOTAL_JOB" + r" +: +([0-9\.]+)" * 4
        re_time = re.compile(expr, re.M)
        with open(filename, "rb") as fmess:
            content = fmess.read().decode(errors="replace")
        for mat in re_time.finditer(content):
            elaps += float(mat.group(3))
    return elaps


# helper functions
def _split(feature):
    parent = feature.split("/")
    base = parent.pop(-1)
    return "/".join(parent), base


def _join(parent, base):
    return parent + "/" + base


def main(argv=None):
    """Execute coverage.

    Arguments:
        argv (list, optional): Command line arguments.

    Returns:
        CoverageAnalysis: object.
    """
    # command arguments parser
    parser = argparse.ArgumentParser(
        usage=USAGE, epilog=EPILOG, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="print more details on stdout")
    parser.add_argument(
        "-f", "--force", action="store_true", help="do not use the current data if exist"
    )
    parser.add_argument("--save", action="store_true", help="save results")
    parser.add_argument(
        "--savetxt",
        action="store_true",
        help="save results as text file (to compare with legacy scripts)",
    )
    parser.add_argument(
        "--limit",
        metavar="N",
        action="store",
        type=int,
        help="report features only covered by a testcase longer than N secs.",
    )
    parser.add_argument(
        "--previous",
        metavar="YYYY-MM-DD",
        action="store",
        default="last",
        help="date of the previous run to compare with",
    )
    parser.add_argument(
        "--wrkdir", type=Path, action="store", default=Path.cwd(), help="working directory"
    )
    parser.add_argument(
        "--srcdir",
        type=Path,
        action="store",
        default=Path.home() / "dev" / "codeaster" / "src",
        help="directory containing code_aster source files",
    )
    parser.add_argument(
        "--installdir",
        type=Path,
        action="store",
        required=True,
        help="installation directory (contains 'bin/run_aster')",
    )
    parser.add_argument(
        "--resdir", type=Path, action="store", help="directory containing testcases results"
    )

    argv = argv or sys.argv[1:]
    args = parser.parse_args(argv)

    os.makedirs(args.wrkdir, exist_ok=True)
    cov = CoverageAnalysis(
        root=args.wrkdir,
        srcdir=args.srcdir,
        installdir=args.installdir,
        resdir=args.resdir,
        previous=args.previous,
        limit=args.limit,
        force=args.force,
        verbose=args.verbose,
    )
    cov.check()
    cov.check_diff()
    if args.save:
        cov.save()
    if args.savetxt:
        cov.savetxt()
    return cov


if __name__ == "__main__":
    main()
