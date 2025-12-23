# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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

from code_aster.Commands import *
from code_aster import CA


# how to compute FORC_NODA for load-displacement plot
def compute_average_forc(RESULTAT, groupsOfCells, NOM_CMP, NOM_PARA):

    # project FORC_NODA on Lagrange space - this is not a force but a stress in practice
    hho = CA.HHO(RESULTAT.getModel())

    _mesh = RESULTAT.getMesh()
    dim = _mesh.getDimension()

    indexes = RESULTAT.getIndexes()
    for idx in indexes:
        forc_noda = RESULTAT.getField("FORC_NODA", idx)
        # project using Face DoFs only since restricted on faces
        forc_nodes = hho.projectOnLagrangeSpace(forc_noda, groupsOfCells)
        # use field HHO_VITE since can not overload
        RESULTAT.setField(forc_nodes, "HHO_VITE", idx)

    # to have a force, we have to integrate on the surface
    dim2typ = {2: "1D", 3: "2D"}
    TAB0 = POST_ELEM(
        INTEGRALE=_F(
            GROUP_MA=groupsOfCells, NOM_CMP=NOM_CMP, NOM_CHAM="HHO_VITE", TYPE_MAILLE=dim2typ[dim]
        ),
        RESULTAT=RESULTAT,
    )

    nbNodes = len(_mesh.getNodesFromCells(groupsOfCells))

    # to have an equivalent of FORC_NODA, we divide by the number of nodes
    coefF = FORMULE(VALE=f"INTE_{NOM_CMP}/nbNodes", NOM_PARA=f"INTE_{NOM_CMP}", nbNodes=nbNodes)

    TAB = CALC_TABLE(TABLE=TAB0, ACTION=_F(OPERATION="OPER", FORMULE=coefF, NOM_PARA=NOM_PARA))

    # IMPR_TABLE(UNITE=6, TABLE=TAB)

    return TAB


def compute_loadVSdisplacement(RESULTAT, groupsOfCells, NOM_CMP):

    # average displacement
    TAB_DEPL = POST_RELEVE_T(
        ACTION=_F(
            OPERATION="MOYENNE_ARITH",
            INTITULE="Ouverture",
            RESULTAT=RESULTAT,
            NOM_CHAM="HHO_DEPL",
            GROUP_MA=groupsOfCells,
            NOM_CMP=NOM_CMP,
        )
    )

    coefO = FORMULE(VALE="MOYENNE", NOM_PARA="MOYENNE")

    TAB1_DEPL = CALC_TABLE(
        TABLE=TAB_DEPL, ACTION=_F(OPERATION="OPER", FORMULE=coefO, NOM_PARA="N" + NOM_CMP)
    )

    # average force
    TAB_FORC = compute_average_forc(RESULTAT, groupsOfCells, NOM_CMP, "F" + NOM_CMP)

    # table load vs displacement
    frdeplac = CALC_TABLE(
        TABLE=TAB1_DEPL,
        ACTION=(
            _F(OPERATION="COMB", TABLE=TAB_FORC, NOM_PARA="INST"),
            _F(OPERATION="EXTR", NOM_PARA=("INST", "N" + NOM_CMP, "F" + NOM_CMP)),
        ),
    )

    return frdeplac
