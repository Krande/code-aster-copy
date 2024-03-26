#! /usr/bin/env python
# encoding: utf-8
# DC 2008
# Thomas Nagy 2016-2018 (ita)
from __future__ import annotations

import re
import pathlib
from dataclasses import dataclass, field

INC_REGEX = r"""(?:^|['">]\s*;)\s*(?:|#\s*)INCLUDE\s+(?:\w+_)?[<"'](.+?)(?=["'>])"""
USE_REGEX = r"""(?:^|;)\s*USE(?:\s+|(?:(?:\s*,\s*(?:NON_)?INTRINSIC)?\s*::))\s*(\w+)"""
MOD_REGEX = r"""(?:^|;)\s*MODULE(?!\s+(?:PROCEDURE|SUBROUTINE|FUNCTION))\s+(\w+)"""
SMD_REGEX = r"""(?:^|;)\s*SUBMODULE\s*\(([\w:]+)\)\s*(\w+)"""

re_inc = re.compile(INC_REGEX, re.I | re.M)
re_use = re.compile(USE_REGEX, re.I | re.M)
re_mod = re.compile(MOD_REGEX, re.I | re.M)
re_smd = re.compile(SMD_REGEX, re.I | re.M)


@dataclass
class fortran_parser:
    incpaths: list[pathlib.Path]
    srcpaths: list[pathlib.Path]
    waiting: list[pathlib] = field(default_factory=list)
    nodes: list[Node] = field(default_factory=list)

    _glob_header_files: list[pathlib.Path] = field(default_factory=list)
    _header_name_map: dict[str, Node] = field(default_factory=dict)
    _processed_files: list[pathlib.Path] = field(default_factory=list)

    def __post_init__(self):
        self._headers = set()
        self._f_files = set()
        self._header_name_map = {}
        for inc_path in self.incpaths:
            if isinstance(inc_path, str):
                inc_path = pathlib.Path(inc_path)
            self._headers.update(set(inc_path.rglob("*.h")))
            self._header_name_map.update({x.relative_to(inc_path).as_posix(): x for x in self._headers})

        for src_path in self.srcpaths:
            if isinstance(src_path, str):
                src_path = pathlib.Path(src_path)
            self._f_files.update(set(src_path.rglob("*.F90")))

        if len(self._headers) != len(self._header_name_map):
            # find duplicates
            names = [x.name for x in self._glob_header_files]
            duplicates = set([x for x in names if names.count(x) > 1])
            raise ValueError(f"Duplicate header files found: {duplicates}")

    def find_deps(self, ffile: pathlib.Path) -> tuple[list[str], list[str], list[str]]:
        """
        Parses a Fortran file to obtain the dependencies used/provided

        :param ffile: fortran file to read
        :type ffile: :py:class:`waflib.Node.Node`
        :return: lists representing the includes, the modules used, and the modules created by a fortran file
        :rtype: tuple of list of strings
        """
        txt = ffile.read_text(encoding='utf-8')

        # Search the entire text for each pattern
        incs = re_inc.findall(txt)
        uses = re_use.findall(txt)
        mods = re_mod.findall(txt)

        # Submodule processing requires special handling to combine the parent module and submodule names
        smd_matches = re_smd.findall(txt)
        for parent, submodule in smd_matches:
            uses.append(parent)  # Assuming you want to track submodule parent modules as a use-case
            mods.append(f"{parent}:{submodule}")

        return incs, uses, mods

    def start(self):
        """
        Start parsing. Use the stack ``self.waiting`` to hold nodes to iterate on
        """
        self.waiting = [x for x in self._f_files]

        while self.waiting:
            nd = self.waiting.pop(0)
            self.iter(nd)

    def iter(self, node):
        """
        Processes a single file during dependency parsing. Extracts files used
        modules used and modules provided.
        """
        incs, uses, mods = self.find_deps(node)
        for x in incs:
            if x in self.seen:
                continue
            self.seen.append(x)
            self.tryfind_header(x)

        for x in uses:
            name = "USE@%s" % x
            if not name in self.names:
                self.names.append(name)

        for x in mods:
            name = "MOD@%s" % x
            if not name in self.names:
                self.names.append(name)

    def tryfind_header(self, filename):
        """
        Adds an include file to the list of nodes to process

        :param filename: file name
        :type filename: string
        """
        found = None
        for n in self.incpaths:
            found = n.find_resource(filename)
            if found:
                self.nodes.append(found)
                self.waiting.append(found)
                break
        if not found:
            if not filename in self.names:
                self.names.append(filename)

    def find_header_resource(self, filename):
        return self._header_name_map.get(filename)


@dataclass
class Node:
    target: pathlib.Path
    includes: list[Node] = field(default_factory=list)
    uses: list[Node] = field(default_factory=list)
    mods: list[Node] = field(default_factory=list)
    smds: list[Node] = field(default_factory=list)

    def read(self):
        return self.target.read_text(encoding="utf8")


if __name__ == "__main__":
    include_dir = pathlib.Path("../bibfor/include").resolve().absolute()
    lib_dir = pathlib.Path("../bibfor").resolve().absolute()

    fparser = fortran_parser([include_dir], [lib_dir])
    fparser.start()
    print(fparser)
