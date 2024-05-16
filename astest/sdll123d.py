# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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

import numpy as NP


def EXTR_MATR(matr, dummy_arg):
    """fonction permettant de récuperer les matrices assemblées au format numpy"""

    matrix = matr.toNumpy()
    dofnum = matr.getDOFNumbering()
    lagr_dofs = dofnum.getLagrangeDOFs()
    phys_dofs = set(dofnum.getPhysicalDOFs())

    # suppression des ddls physiques
    remove_dofs = set()
    for lagr_dof in lagr_dofs:
        node, comp = dofnum.getNodeAndComponentFromDOF(lagr_dof)
        comp = comp.split(":")[1]
        for phys_dof in phys_dofs:
            if dofnum.getNodeAndComponentFromDOF(phys_dof) == (node, comp):
                remove_dofs.add(phys_dof)
    phys_dofs -= remove_dofs
    phys_dofs = NP.array(list(phys_dofs))
    matrix = matrix[phys_dofs, :][:, phys_dofs]

    # noms des ddls physiques
    dof_names = []
    for phys_dof in phys_dofs:
        dof_names.append(dofnum.getComponentFromDOF(phys_dof))

    return dof_names, matrix
