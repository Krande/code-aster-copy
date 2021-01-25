# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

# person_in_charge: mathieu.courtois@edf.fr
"""
:py:class:`FieldOnNodesReal` --- Fields defined on nodes of elements
**********************************************************************
"""

import numpy

import aster
from libaster import FieldOnNodesReal, DOFNumbering

from ..Utilities import injector


@injector(FieldOnNodesReal)
class ExtendedFieldOnNodesReal(object):
    cata_sdj = "SD.sd_champ.sd_cham_no_class"

    def EXTR_COMP(self, comp=' ', lgno=[], topo=0):
        """ retourne les valeurs de la composante comp du champ sur la liste
        de groupes de noeuds lgno avec eventuellement l'info de la
        topologie si topo>0. Si lgno est une liste vide, c'est equivalent
        a un TOUT='OUI' dans les commandes aster

        Arguments:
            comp (str): Name of the component.
            lgno (list[str]): List of groups of nodes.
            topo (int): ``topo == 1`` means to return informations about the
                support of the field.

        Returns:
            :py:class:`.post_comp_cham_no`: Object containing the values and,
            eventually, the topological informations of the support.
        """

        ncham = self.getName()
        ncham = ncham + (19 - len(ncham)) * ' '
        nchams = aster.get_nom_concept_unique('_')
        ncmp = comp + (8 - len(comp)) * ' '

        aster.prepcompcham(ncham, nchams, ncmp, "NO      ", topo, lgno)

        valeurs = numpy.array(
            aster.getvectjev(nchams + (19 - len(nchams)) * ' ' + '.V'))

        if (topo > 0):
            noeud = (
                aster.getvectjev(nchams + (19 - len(nchams)) * ' ' + '.N'))
        else:
            noeud = None

        if comp[:1] == ' ':
            comp = (aster.getvectjev(nchams + (19 - len(nchams)) * ' ' + '.C'))
            aster.prepcompcham(
                "__DETR__", nchams, ncmp, "NO      ", topo, lgno)
            return post_comp_cham_no(valeurs, noeud, comp)
        else:
            aster.prepcompcham(
                "__DETR__", nchams, ncmp, "NO      ", topo, lgno)
            return post_comp_cham_no(valeurs, noeud)

    def setDirichletBC(self, **kwargs):
        """Set the values of the Dirichlet boundary conditions of the degrees
        of freedom on the group of nodes or cells.

        Arguments:
            GROUP_MA (str or list(str)): name of the group of cells
            GROUP_NO (str or list(str)): name of the group of nodes.
            NOEUD    (int or list(int)): indices of the nodes.
            DX, DY, ... (float): name and value of the degree of freedom
        """
        ldofNumbering = [dep for dep in self.getDependencies() if isinstance(dep,DOFNumbering)]
        if not ldofNumbering:
            raise RuntimeError("Cannot retrieve dofNumbering")
        dofNumbering = [dep for dep in self.getDependencies() if isinstance(dep,DOFNumbering)][-1]
        mesh = dofNumbering.getMesh()
        if mesh.isParallel():
            raise RuntimeError("No support for ParallelDOFNumbering")
        # build the indirection table between (nodeid, dof) and row
        indir = {}
        for row in dofNumbering.getRowsAssociatedToLagrangeMultipliers():
            dof = dofNumbering.getComponentAssociatedToRow(row)
            node = -1 * int(dofNumbering.getNodeAssociatedToRow(row))  # constrained nodes have id < 0
            if node > 0 and dof != "":
                indir.setdefault((node,dof), []).append(row)  # there may be 2 Lagrange multipliers per constraint
        # build the group of nodes to be processed
        if not ('GROUP_MA' in kwargs.keys() or
                'GROUP_NO' in kwargs.keys() or
                'NOEUD' in kwargs.keys()):
            raise ValueError("a group of cells (GROUP_MA), a group of nodes (GROUP_NO)"
                             " or a list of nodes (NOEUD) must be provided")
        lNodes = []
        if 'GROUP_MA' in kwargs.keys():
            lGrpMa = kwargs['GROUP_MA']
            lGrpMa = [lGrpMa] if isinstance(lGrpMa, str) else lGrpMa
            for grMa in lGrpMa:
                if mesh.hasGroupOfNodes(grMa):
                    nodes = mesh.getNodes(grMa)
                    lNodes += nodes
                elif mesh.hasGroupOfCells(grMa):
                    connec = mesh.getConnectivity()
                    nodes = [node for cell in mesh.getCells(grMa) for node in connec[cell - 1]]
                    lNodes += nodes
                else:
                    raise ValueError("no {} group of cells".format(grMa))
        if 'GROUP_NO' in kwargs.keys():
            lgrNo = kwargs['GROUP_NO']
            lgrNo = [lgrNo] if isinstance(lgrNo, str) else lgrNo
            for grNo in lgrNo:
                if mesh.hasGroupOfNodes(grNo):
                    nodes = mesh.getNodes(grNo)
                    lNodes += nodes
                else:
                    raise ValueError("no {} group of nodes".format(grNo))
        if 'NOEUD' in kwargs.keys():
            lNo = kwargs['NOEUD']
            lNo = [lNo] if isinstance(lNo, int) else lNo
            lNodes += lNo
        # keep unique node id
        lNodes = list(set(lNodes))
        # change the given values
        assignedDOF=0
        self.updateValuePointers()  # update Jeveux pointers before assignment
        for node in lNodes:
            for (dof, val) in kwargs.items():
                if dof in ['GROUP_MA','GROUP_NO','NOEUD']:  # only process DOF here
                    continue
                if (node, dof) in indir.keys():
                    assignedDOF += 1
                    for row in indir[(node, dof)]:
                        self[row-1] = val   # row are 1-based
        if assignedDOF == 0:
            raise ValueError("No bounday condition has been set - no entity handle the given degree of freedom")




class post_comp_cham_no:
    """Container object that store the results of
    :py:meth:`code_aster.Objects.FieldOnNodesReal.EXTR_COMP`.

    The support of the field may be unknown. In this case, :py:attr:`noeud`
    and :py:attr:`comp` are set to *None*.

    Attributes:
        valeurs (numpy.ndarray): Values of the field.
        noeud (list[int]): List of nodes numbers.
        comp (list[int]): List of components.
    """

    def __init__(self, valeurs, noeud=None, comp=None):
        self.valeurs = valeurs
        self.noeud = noeud
        self.comp = tuple(i.strip() for i in comp) if comp else None
