# coding: utf-8

# Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
#
# This file is part of Code_Aster.
#
# Code_Aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Code_Aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.


from ..Helpers import FileAccess, LogicalUnitFile
from ..Messages import UTMESS
from ..Objects import Mesh, ParallelMesh
from ..Supervis import ExecuteCommand
from ..Utilities import haveMPI

from .pre_gibi import PRE_GIBI
from .pre_gmsh import PRE_GMSH
from .pre_ideas import PRE_IDEAS


class MeshReader(ExecuteCommand):
    """Command that creates a :class:`~code_aster.Objects.Mesh` from a file."""

    command_name = "LIRE_MAILLAGE"

    def _create_parallel_mesh(self, keywords):
        """Tell if the command is creating a ParallelMesh

        Arguments:
            keywords (dict): User's keywords.

        Returns:
            bool: *True* if a ParallelMesh is creating, *False* otherwise.
        """
        return keywords["FORMAT"] == "MED" and keywords["PARTITIONNEUR"] == "PTSCOTCH" and haveMPI()

    def create_result(self, keywords):
        """Create the :class:`~code_aster.Objects.Mesh`.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        if self._create_parallel_mesh(keywords):
            self._result = ParallelMesh()
        else:
            self._result = Mesh()

    def exec_(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """
        need_conversion = ("GIBI", "GMSH", "IDEAS")
        unit = keywords.get("UNITE")
        fmt = keywords["FORMAT"]

        if fmt in need_conversion:
            tmpfile = LogicalUnitFile.new_free(access=FileAccess.New)
            unit_op = tmpfile.unit
            keywords["FORMAT"] = "ASTER"
        else:
            unit_op = unit

        if fmt == "GIBI":
            PRE_GIBI(UNITE_GIBI=unit, UNITE_MAILLAGE=unit_op)
        elif fmt == "GMSH":
            PRE_GMSH(UNITE_GMSH=unit, UNITE_MAILLAGE=unit_op)
        elif fmt == "IDEAS":
            coul = keywords.pop("CREA_GROUP_COUL", "NON")
            PRE_IDEAS(UNITE_IDEAS=unit, UNITE_MAILLAGE=unit_op, CREA_GROUP_COUL=coul)

        assert keywords["FORMAT"] in ("ASTER", "MED")
        if keywords["FORMAT"] == "ASTER":
            if fmt in need_conversion:
                tmpfile.release()
            keywords["UNITE"] = unit_op
            super().exec_(keywords)
            return

        meshname = keywords.get("NOM_MED")
        verbose = keywords["INFO"]
        if keywords.get("INFO_MED", 0):
            # cheat code for debugging and detailed time informations
            verbose |= 4
        if meshname is None:
            meshname = ""
        if self._result.isParallel():
            filename = LogicalUnitFile.filename_from_unit(unit)
            self._result.readMedFile(filename, meshname, partitioned=False, verbose=verbose)
        else:
            if keywords["PARTITIONNEUR"] == "PTSCOTCH":
                assert not haveMPI()
                UTMESS("A", "MED_18")
            filename = LogicalUnitFile.filename_from_unit(unit)
            self._result.readMedFile(filename, meshname, verbose=verbose)

        if keywords["VERI_MAIL"]["VERIF"] == "OUI":
            self._result.check(keywords["VERI_MAIL"]["APLAT"])

    def post_exec(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """
        self._result.build()


LIRE_MAILLAGE = MeshReader.run
