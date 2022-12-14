# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
********************************************************************

The *Field On Nodes* object exists for real numbers (:py:class:`FieldOnNodesReal`),
integers (:py:class:`FieldOnNodesLong`),
strings (:py:class:`FieldOnNodesChar8`) and
complex numbers (:py:class:`FieldOnNodesComplex`).
"""

import numpy, functools

import aster
from libaster import (
    FieldOnNodesReal,
    FieldOnNodesLong,
    FieldOnNodesChar8,
    FieldOnNodesComplex,
    DOFNumbering,
)
from ..Objects.Serialization import InternalStateBuilder
from ..Utilities import injector


class FieldOnNodesStateBuilder(InternalStateBuilder):
    """Class that returns the internal state of a *FieldOnNodes*."""

    def save(self, field):
        """Return the internal state of a *Result* to be pickled.

        Arguments:
            field (*FieldOnNodes*): The *FieldOnNodes* object to be pickled.

        Returns:
            *InternalStateBuilder*: The internal state itself.
        """
        super().save(field)
        self._st["mesh"] = field.getMesh()
        self._st["dofd"] = field.getDescription()
        return self

    def restore(self, field):
        """Restore the *DataStructure* content from the previously saved internal
        state.

        Arguments:
            field (*DataStructure*): The *DataStructure* object to be pickled.
        """
        super().restore(field)
        if self._st["mesh"]:
            field.setMesh(self._st["mesh"])
        if self._st["dofd"]:
            field.setDescription(self._st["dofd"])


@injector(FieldOnNodesReal)
class ExtendedFieldOnNodesReal:
    cata_sdj = "SD.sd_champ.sd_cham_no_class"
    internalStateBuilder = FieldOnNodesStateBuilder

    def EXTR_COMP(self, comp=" ", lgno=[], topo=0):
        """retourne les valeurs de la composante comp du champ sur la liste
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
        ncham = ncham + (19 - len(ncham)) * " "
        nchams = aster.get_nom_concept_unique("_")
        ncmp = comp + (8 - len(comp)) * " "

        aster.prepcompcham(ncham, nchams, ncmp, "NO      ", topo, lgno)

        valeurs = numpy.array(aster.getvectjev(nchams + (19 - len(nchams)) * " " + ".V"))

        if topo > 0:
            noeud = aster.getvectjev(nchams + (19 - len(nchams)) * " " + ".N")
        else:
            noeud = None

        if comp[:1] == " ":
            comp = aster.getvectjev(nchams + (19 - len(nchams)) * " " + ".C")
            aster.prepcompcham("__DETR__", nchams, ncmp, "NO      ", topo, lgno)
            return post_comp_cham_no(valeurs, noeud, comp)
        else:
            aster.prepcompcham("__DETR__", nchams, ncmp, "NO      ", topo, lgno)
            return post_comp_cham_no(valeurs, noeud)

    @property
    @functools.lru_cache()
    def __NodeDOF2Row(self):
        """Build the indirection table between (nodeid, dof) and row"""
        ldofNumbering = [dep for dep in self.getDependencies() if isinstance(dep, DOFNumbering)]
        if not ldofNumbering:
            raise RuntimeError("Cannot retrieve dofNumbering")
        dofNumbering = [dep for dep in self.getDependencies() if isinstance(dep, DOFNumbering)][-1]
        # build the indirection table between (nodeid, dof) and row
        indir = {}
        for row in dofNumbering.getRowsAssociatedToLagrangeMultipliers():
            if not dofNumbering.isRowAssociatedToPhysical(row):
                dof = dofNumbering.getComponentAssociatedToRow(row).split(":")[-1]
                if dof != "MPC":
                    node = dofNumbering.getNodeAssociatedToRow(row)
                    # there may be 2 Lagrange multipliers per constraint
                    indir.setdefault((node, dof), []).append(row)
        return indir

    def setDirichletBC(self, **kwargs):
        """Set the values of the Dirichlet boundary conditions of the degrees
        of freedom on the group of nodes or cells.

        Arguments:
            GROUP_MA (str or list(str)): name of the group of cells
            GROUP_NO (str or list(str)): name of the group of nodes.
            NOEUD    (int or list(int)): indices of the nodes.
            DX, DY, ... (float): name and value of the degree of freedom
        """
        ldofNumbering = [dep for dep in self.getDependencies() if isinstance(dep, DOFNumbering)]
        if not ldofNumbering:
            raise RuntimeError("Cannot retrieve dofNumbering")
        dofNumbering = [dep for dep in self.getDependencies() if isinstance(dep, DOFNumbering)][-1]
        mesh = dofNumbering.getMesh()
        if mesh.isParallel():
            raise RuntimeError("No support for ParallelDOFNumbering")
        # build the group of nodes to be processed
        if not (
            "GROUP_MA" in kwargs.keys() or "GROUP_NO" in kwargs.keys() or "NOEUD" in kwargs.keys()
        ):
            raise ValueError(
                "a group of cells (GROUP_MA), a group of nodes (GROUP_NO)"
                " or a list of nodes (NOEUD) must be provided"
            )
        lNodes = []
        if "GROUP_MA" in kwargs.keys():
            lGrpMa = kwargs["GROUP_MA"]
            lGrpMa = [lGrpMa] if isinstance(lGrpMa, str) else lGrpMa
            connec = mesh.getConnectivity()
            for grMa in lGrpMa:
                if mesh.hasGroupOfNodes(grMa):
                    nodes = mesh.getNodes(grMa)
                    lNodes += nodes
                elif mesh.hasGroupOfCells(grMa):
                    nodes = [node-1 for cell in mesh.getCells(grMa) for node in connec[cell]]
                    lNodes += nodes
                else:
                    raise ValueError("no {} group of cells".format(grMa))
        if "GROUP_NO" in kwargs.keys():
            lgrNo = kwargs["GROUP_NO"]
            lgrNo = [lgrNo] if isinstance(lgrNo, str) else lgrNo
            for grNo in lgrNo:
                if mesh.hasGroupOfNodes(grNo):
                    nodes = mesh.getNodes(grNo)
                    lNodes += nodes
                else:
                    raise ValueError("no {} group of nodes".format(grNo))
        if "NOEUD" in kwargs.keys():
            lNo = kwargs["NOEUD"]
            lNo = [lNo-1] if isinstance(lNo, int) else [no-1 for no in lNo]
            lNodes += lNo
        # keep unique node id
        lNodes = list(set(lNodes))
        # change the given values
        assignedDOF = 0
        self.updateValuePointers()  # update Jeveux pointers before assignment
        for node in lNodes:
            for (dof, val) in kwargs.items():
                if dof in ["GROUP_MA", "GROUP_NO", "NOEUD"]:  # only process DOF here
                    continue
                if (node, dof) in self.__NodeDOF2Row.keys():
                    assignedDOF += 1
                    for row in self.__NodeDOF2Row[(node, dof)]:
                        self[row] = val
        if assignedDOF == 0:
            raise ValueError(
                "No bounday condition has been set - no entity handle the given degree of freedom"
            )


@injector(FieldOnNodesLong)
class ExtendedFieldOnNodesLong:
    cata_sdj = "SD.sd_champ.sd_cham_no_class"
    internalStateBuilder = FieldOnNodesStateBuilder


@injector(FieldOnNodesChar8)
class ExtendedFieldOnNodesChar8:
    cata_sdj = "SD.sd_champ.sd_cham_no_class"
    internalStateBuilder = FieldOnNodesStateBuilder


@injector(FieldOnNodesComplex)
class ExtendedFieldOnNodesComplex:
    cata_sdj = "SD.sd_champ.sd_cham_no_class"
    internalStateBuilder = FieldOnNodesStateBuilder

    def EXTR_COMP(self, comp=" ", lgno=[], topo=0):
        """retourne les valeurs de la composante comp du champ sur la liste
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
        ncham = ncham + (19 - len(ncham)) * " "
        nchams = aster.get_nom_concept_unique("_")
        ncmp = comp + (8 - len(comp)) * " "

        aster.prepcompcham(ncham, nchams, ncmp, "NO      ", topo, lgno)

        valeurs = numpy.array(aster.getvectjev(nchams + (19 - len(nchams)) * " " + ".V"))

        if topo > 0:
            noeud = aster.getvectjev(nchams + (19 - len(nchams)) * " " + ".N")
        else:
            noeud = None

        if comp[:1] == " ":
            comp = aster.getvectjev(nchams + (19 - len(nchams)) * " " + ".C")
            aster.prepcompcham("__DETR__", nchams, ncmp, "NO      ", topo, lgno)
            return post_comp_cham_no(valeurs, noeud, comp)
        else:
            aster.prepcompcham("__DETR__", nchams, ncmp, "NO      ", topo, lgno)
            return post_comp_cham_no(valeurs, noeud)


class post_comp_cham_no:
    """Container object that store the results of
    :py:meth:`code_aster.Objects.FieldOnNodesReal.EXTR_COMP` and
    :py:meth:`code_aster.Objects.FieldOnNodesComplex.EXTR_COMP`.

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
