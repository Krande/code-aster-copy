import re
import pathlib

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent

INC_REGEX = r"""(?:^|['">]\s*;)\s*(?:|#\s*)INCLUDE\s+(?:\w+_)?[<"'](.+?)(?=["'>])"""
USE_REGEX = r"""(?:^|;)\s*USE(?:\s+|(?:(?:\s*,\s*(?:NON_)?INTRINSIC)?\s*::))\s*(\w+)"""
MOD_REGEX = r"""(?:^|;)\s*MODULE(?!\s+(?:PROCEDURE|SUBROUTINE|FUNCTION))\s+(\w+)"""
SMD_REGEX = r"""(?:^|;)\s*SUBMODULE\s*\(([\w:]+)\)\s*(\w+)"""

re_inc = re.compile(INC_REGEX, re.I)
re_use = re.compile(USE_REGEX, re.I)
re_mod = re.compile(MOD_REGEX, re.I)
re_smd = re.compile(SMD_REGEX, re.I)

for_include_dir = ROOT_DIR / "bibfor/include"
c_include_dir = ROOT_DIR / "bibc/include"
cxx_include_dir = ROOT_DIR / "bibcxx/include"


def resolve_dependency_order(ffile1: pathlib.Path):
    if not ffile1.exists():
        print(f"File {file1} does not exist")
    txt = ffile1.read_text()
    incs = []
    uses = []
    mods = []
    for line in txt.splitlines():
        # line by line regexp search? optimize?
        m = re_inc.search(line)
        if m:
            incs.append(m.group(1))
        m = re_use.search(line)
        if m:
            uses.append(m.group(1))
        m = re_mod.search(line)
        if m:
            mods.append(m.group(1))
        m = re_smd.search(line)
        if m:
            mods.append(m.group(2))
    return incs, uses, mods


if __name__ == '__main__':
    file1 = ROOT_DIR / "bibfor/utilitai/utgtme.F90"

    result = resolve_dependency_order(file1)
    print(result)
