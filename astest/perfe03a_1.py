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
from code_aster.Commands import DEFI_GROUP
from random import random
import numpy as NP


def F_DEFI_GROUP(MAIL):
    def dist(x, y):
        """computes the square distance between 2 points"""
        v = y - x
        return NP.sqrt(NP.dot(v, v))

    ### procedure of definition of the sets of grains within the mesh
    def DEFI_GRAINS(nb_germs, dim, prefix_nomgr, mesh_name="MAIL", **args):
        """We here introduce a certain number of germs for grains. Each
        germ grows until it reaches the area delimited by other germs. All
        elements sets (GROUP_MA) belonging to the same germ are identified
        and returned. The germ to which belong a element set is evaluated
        as the one whose distance to middle of the elset is minimum.
        - nb_germs corresponds to the number of germs
        - mesh_name corresponds to the name of the Aster mesh
        - LISTE_GROUP is a keyword CREA_GROUP_MA of DEFI_GROUP produced as
        output, which contains the definition of the grains.
        - dim is the dimension of the cubic r.e.v.
        """

        # definition of germs
        germ0 = NP.array([dim * random(), dim * random(), dim * random()])
        germs = [germ0]

        #  computation of the 'repulsion distance'
        V = dim**3
        rep_dist = ((V / nb_germs) ** (1.0 / 3.0)) * 0.8
        print("Number of germs: %i" % (nb_germs,))
        print("Minimal distance between 2 germs: %5.4e" % (rep_dist,))

        # computation of the random coordinates of all germs with respect to the
        # repulsion distance
        print("Computing the coordinates of the germs...")
        for i in range(1, nb_germs):
            good_point = False
            while not good_point:
                this_germ = NP.array([dim * random(), dim * random(), dim * random()])
                good_point = True
                for germ in germs:
                    if NP.sqrt(NP.dot((germ - this_germ), (germ - this_germ))) < rep_dist:
                        good_point = False
            germs.append(this_germ)
        print("Done.")

        # get the coordinates of the mesh nodes, and the number of nodes
        print("Getting mesh nodes and coordinates...")
        full_mesh_name = "%-8s" % (mesh_name,)
        coordo = MAIL.getCoordinates().getValues()
        nb_nodes = len(coordo) / 3
        coordo = NP.reshape(coordo, [nb_nodes, 3])
        print("Done.")

        # get the mesh connectivity, the elements sets and names, the node names
        print("Getting mesh connectivity and element sets...")
        connectivity = MAIL.getConnextivity()
        type_3D = [18, 19, 20, 21, 22, 23, 24, 25, 26]
        # mailles 3D :
        # 18 : TETRA3
        # 19 : TETRA10
        # 20 : PENTA6',
        # 21 : PENTA15',
        # 22 : PYRAM5',
        # 23 : PYRAM13','
        # 24 : HEXA8',  '
        # 25 : HEXA20', '
        # 26 : HEXA27'
        nb_elset = MAIL.getNumberOfCells()
        listegroup = [None] * nb_elset
        print("Done.")

        # building of the listegroup table indicating to which germ belong each
        # elset
        for i_elset in range(nb_elset):
            ### Warning: in code_aster connectivity table, nodes number begin at 1
            ### whereas in python tables they begin at 0!
            this_elset_nodes = connectivity[i_elset]
            # 3D element only
            if MAIL.getCellType(i_elset) in type_3D:
                ### computation of the barycentre of the elset
                middle_point = NP.array([0.0, 0.0, 0.0])
                nbnodel = len(this_elset_nodes)
                for i in range(nbnodel):
                    middle_point += (1.0 / nbnodel) * coordo[this_elset_nodes[i] - 1]
                ### computation of the minimum distance of the middle point
                ### to all germs
                min_dist = 1.0e10
                index_min = 0
                for i in range(nb_germs):
                    this_germ = germs[i]
                    this_dist = NP.sqrt(
                        NP.dot((this_germ - middle_point), (this_germ - middle_point))
                    )
                    if this_dist < min_dist:
                        min_dist, index_min = this_dist, i
                listegroup[i_elset] = index_min

        # we create now a list of factor keyword to be used in DEFI_GROUP
        # each element is a dictionnary containing the name of the grain
        # (keyword 'NOM') and a list of elsets names (keyword 'MAILLE')
        LIST_GROUP = []
        for i_germ in range(nb_germs):
            this_dict = {}
            this_dict["NOM"] = prefix_nomgr + str(i_germ + 1)
            this_elset_list = []
            for i_elset in range(nb_elset):
                if listegroup[i_elset] == i_germ:
                    this_elset_list.append(MAIL.getNodeName(i_elset))
            this_dict["MAILLE"] = this_elset_list
            LIST_GROUP.append(this_dict)

        # safely return LIST_GROUP
        return LIST_GROUP

    print("Definition of grains...")

    LISTE_GROUP = DEFI_GRAINS(Ng, dim, prefix_nomgr)

    print("Grains done")
    print("Adding the grains to the mesh...")

    MAIL = DEFI_GROUP(reuse=MAIL, MAILLAGE=MAIL, CREA_GROUP_MA=LISTE_GROUP)

    print("Grains added")

    return MAIL
