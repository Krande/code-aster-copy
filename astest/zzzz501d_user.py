# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------


class SimpleUserObject:
    """'Hello World!' of embedding code_aster objects."""

    aster_embedded = ["attrname"]

    def __init__(self, obj):
        self.attrname = obj

    def __getstate__(self):
        return [self.attrname]

    def __setstate__(self, state):
        assert len(state) == 1, state
        self.attrname = state[0]


class ComplexUserObject:
    """Example for user class that embeds code_aster objects in its attributes."""

    def __init__(self, mesh, resu, values):
        self.asattr = mesh
        self.aslist = [mesh, resu]
        self.astuple = (mesh, resu)
        self.asdict = {"mesh": mesh, "result": resu}
        self.nested = {"lotod": ["list", ("tuple", {"dict": mesh})]}
        self.subobj = SimpleUserObject(mesh)
        self.values = values

    def __getstate__(self):
        return [
            self.asattr,
            self.aslist,
            self.astuple,
            self.asdict,
            self.nested,
            self.subobj,
            self.values,
        ]

    def __setstate__(self, state):
        assert len(state) == 7, state
        (
            self.asattr,
            self.aslist,
            self.astuple,
            self.asdict,
            self.nested,
            self.subobj,
            self.values,
        ) = state

    def __repr__(self):
        """Representation for debugging."""
        lines = [
            ">>> Object content <<<",
            f"asattr: {self.asattr}",
            f"aslist: {self.aslist}",
            f"astuple: {self.astuple}",
            f"asdict: {self.asdict}",
            f"nested: {self.nested}",
            f"subobj: {self.subobj}",
            f"values: {self.values}",
        ]
        return "\n".join(lines)
