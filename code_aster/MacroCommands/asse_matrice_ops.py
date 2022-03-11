# coding: utf-8

# Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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

from libaster import setFortranLoggingLevel, resetFortranLoggingLevel

from ..Objects import (
    AssemblyMatrixDisplacementReal,
    AssemblyMatrixDisplacementComplex,
    AssemblyMatrixTemperatureReal,
    AssemblyMatrixPressureComplex,
    ElementaryMatrixDisplacementReal,
    ElementaryMatrixDisplacementComplex,
    ElementaryMatrixTemperatureReal,
    ElementaryMatrixPressureComplex
)
from ..Utilities import force_list


def asse_matrice_ops(self, **args):
    """Execute the command ASSE_MATRICE.

    Arguments:
        **args (dict): User's keywords.

    Returns:
        AssemblyMatrix: matrix assemblyed from MATR_ELEM
    """

    verbosity = args["INFO"]
    setFortranLoggingLevel(verbosity)

    # Create result
    matr_elem = force_list(args["MATR_ELEM"])
    if (isinstance(matr_elem[0], ElementaryMatrixDisplacementReal)):
        matr = AssemblyMatrixDisplacementReal()
    elif(isinstance(matr_elem[0], ElementaryMatrixDisplacementComplex)):
        matr = AssemblyMatrixDisplacementComplex()
    elif(isinstance(matr_elem[0], ElementaryMatrixTemperatureReal)):
        matr = AssemblyMatrixTemperatureReal()
    elif(isinstance(matr_elem[0], ElementaryMatrixPressureComplex)):
        matr = AssemblyMatrixPressureComplex()
    else:
        raise RuntimeError("Unexpected matrix")

    # Add MATR_ELEM
    for elem in matr_elem:
        if isinstance(elem, type(matr_elem[0])):
            matr.addElementaryMatrix(elem)
        else:
            raise RuntimeError("Incompatible elementary matrices")

    # Add NUME_DDL
    matr.setDOFNumbering(args["NUME_DDL"])

    # Add DirichletBC
    if "CHAR_CINE" in args:
        for cine in args["CHAR_CINE"]:
            matr.addDirichletBC(cine)

    # Assemble
    matr.assemble(False)

    if args.get("SYME") == "OUI":
        matr.symmetrize()

    resetFortranLoggingLevel()

    return matr
