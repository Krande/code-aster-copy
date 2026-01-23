# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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


def barycenter(coor, nodes):

    bar = [0.0, 0.0, 0.0]
    for n in nodes:
        coor_n = coor.getNode(n)
        bar[0] += coor_n.x()
        bar[1] += coor_n.y()
        bar[2] += coor_n.z()

    bar[0] /= float(len(nodes))
    bar[1] /= float(len(nodes))
    bar[2] /= float(len(nodes))

    return bar


def node_up(coor, connec, cell_id):
    bar = barycenter(coor, connec[cell_id])

    if bar[1] > 0.0:
        return True

    return False


def create_grp_2d(mesh):

    connec = mesh.getConnectivity()
    coor = mesh.getCoordinates()

    all_cells = mesh.getCells("Sym")
    nodes_haut = []
    nodes_bas = []

    for c_id in all_cells:
        if node_up(coor, connec, c_id):
            nodes_haut.append(f"M{c_id+1}")
        else:
            nodes_bas.append(f"M{c_id+1}")

    mesh = DEFI_GROUP(
        MAILLAGE=mesh,
        reuse=mesh,
        CREA_GROUP_MA=(
            _F(NOM="LIAISON_Haut", GROUP_MA=("Group_2")),
            _F(NOM="LIAISON_Bas", GROUP_MA=("Group_1")),
            _F(NOM="Encast", GROUP_MA=("Group_3")),
            _F(NOM="Pres", GROUP_MA=("Group_4")),
            _F(NOM="Sym_Haut", MAILLE=nodes_haut),
            _F(NOM="Sym_Bas", MAILLE=nodes_bas),
        ),
    )

    return mesh
