from __future__ import annotations

import pathlib
import re
from dataclasses import dataclass

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent

INC_REGEX = r"""(?:^|['">]\s*;)\s*(?:|#\s*)INCLUDE\s+(?:\w+_)?[<"'](.+?)(?=["'>])"""
USE_REGEX = r"""(?:^|;)\s*USE(?:\s+|(?:(?:\s*,\s*(?:NON_)?INTRINSIC)?\s*::))\s*(\w+)"""
MOD_REGEX = r"""(?:^|;)\s*MODULE(?!\s+(?:PROCEDURE|SUBROUTINE|FUNCTION))\s+(\w+)"""
SMD_REGEX = r"""(?:^|;)\s*SUBMODULE\s*\(([\w:]+)\)\s*(\w+)"""

re_inc = re.compile(INC_REGEX, re.I)
re_use = re.compile(USE_REGEX, re.I)
re_mod = re.compile(MOD_REGEX, re.I)
re_smd = re.compile(SMD_REGEX, re.I)

for_src_dir = ROOT_DIR / "bibfor"
for_include_dir = for_src_dir / "include"

c_src_dir = ROOT_DIR / "bibc"
c_include_dir = c_src_dir / "include"

cxx_src_dir = ROOT_DIR / "bibcxx"
cxx_include_dir = cxx_src_dir / "include"


def resolve_dependencies(ffile1: pathlib.Path) -> list[Node]:
    if not ffile1.exists():
        print(f"File {file1} does not exist")
    txt = ffile1.read_text()
    incs = []

    for line in txt.splitlines():
        # line by line regexp search? optimize?
        m = re_inc.search(line)
        if m:
            incl = m.group(1)
            fortran_include = for_include_dir / incl
            if fortran_include.exists():
                src_file = (for_src_dir / incl).with_suffix('.F90')
                source_files = []
                if src_file.exists():
                    source_files.append(Node(target=src_file, header_dependencies=[], source_dependencies=[]))
                incs.append(Node(target=fortran_include, header_dependencies=[], source_dependencies=source_files))
            elif c_include_dir.exists():
                incs.append(Node(target=c_include_dir / incl, header_dependencies=[], source_dependencies=[]))
            elif cxx_include_dir.exists():
                incs.append(Node(target=cxx_include_dir / incl, header_dependencies=[], source_dependencies=[]))

    return incs


@dataclass
class Node:
    target: pathlib.Path
    header_dependencies: list[Node]
    source_dependencies: list[Node]


def main(start_files: list[pathlib.Path]) -> list[Node]:
    nodes = start_files
    processed_nodes = []
    while nodes:
        node = nodes.pop()
        incs, uses, mods = resolve_dependencies(node)
        for inc in incs:
            if inc not in nodes:
                nodes.append(inc)

        processed_nodes.append(Node(target=node, header_dependencies=incs, source_dependencies=[]))

    return processed_nodes


if __name__ == '__main__':
    file1 = ROOT_DIR / "bibfor/utilitai/utgtme.F90"

    result = main([file1])
    print(result)
