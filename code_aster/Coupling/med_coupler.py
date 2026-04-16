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

"""
Definition of a convenient object to synchronize MEDCoupling fields.
"""

import time

from ..Objects import LoadResult, SimpleFieldOnCellsReal, SimpleFieldOnNodesReal, ParallelMesh
from ..Utilities import ParaMEDMEM as PMM
from ..Utilities import logger, no_new_attributes, MPI
from ..Utilities import medcoupling as MEDC

# need mecoupling >= 9.16.0 to use InterpKernelDECWithOverlap
# remove PMM.InterpKernelDEC later
IKDEC = getattr(PMM, "InterpKernelDECWithOverlap", PMM.InterpKernelDEC)


class CoupledField(PMM.ParaFIELD):
    """Define the properties of an coupled field.

    Attributes:
        support (TypeOfField): Support: ON_NODES or ON_NODES_FE or ON_CELLS.
        td (TypeOfTimeDiscretization): time discretization
        dec (ExtendedInterpKernelDECWithOverlap): dec.
        topo (ComponentTopology): topology.
    """

    def __init__(self, sup, td, dec, topo):
        assert sup in (MEDC.ON_NODES, MEDC.ON_NODES_FE, MEDC.ON_CELLS), sup
        self.dec = dec
        super().__init__(sup, td, dec.mesh, topo)

    def fillWithZero(self):
        self.getArray().fillWithZero()

    def setArray(self, array):
        assert self.getNbOfElems() == array.getNbOfElems()
        self.getField().setArray(array)

    def getNbOfElems(self):
        return self.getArray().getNbOfElems()

    def getInfoOnComponents(self):
        return self.getArray().getInfoOnComponents()

    def setInfoOnComponents(self, info):
        return self.getArray().setInfoOnComponents(info)

    def setDescription(self, desc):
        self.getField().setDescription(desc)

    def getArray(self):
        return self.getField().getArray()

    def getValues(self):
        return self.getArray().getValues()

    def setNature(self, nature):
        return self.getField().setNature(nature)

    def getTypeOfField(self):
        return self.getField().getTypeOfField()

    @property
    def sup(self):
        return self.getField().getTypeOfField()

    @property
    def sup_name(self):
        type_field = self.getField().getTypeOfField()
        if type_field == MEDC.ON_NODES:
            return "nodes"
        elif type_field == MEDC.ON_NODES_FE:
            return "nodes_fe"
        elif type_field == MEDC.ON_CELLS:
            return "cells"
        return "None"


class ExtendedInterpKernelDECWithOverlap(IKDEC):
    """Object that represents a DEC and the necessary properties.

    Arguments:
        src_ranks (list[int]): source procs IDs.
        trg_ranks (list[int]): target procs IDs.
    """

    mesh = _synced = None

    def __init__(self, src_ranks, trg_ranks):
        self.mesh = None
        self._synced = False
        super().__init__(src_ranks, trg_ranks)

    @property
    def synced(self):
        """bool: Tell if the DEC has already been synced."""
        return self._synced

    def synchronize(self):
        """Wrapper on DEC function."""
        if self._synced:
            return
        self._synced = True
        return super().synchronize()


class ExtendedCFEMDEC(PMM.CFEMDEC):
    """Object that represents a DEC and the necessary properties.

    Arguments:
        src_ranks (list[int]): source procs IDs.
        trg_ranks (list[int]): target procs IDs.
    """

    mesh = _synced = None
    pfield = None

    def __init__(self, src_ranks, trg_ranks):
        self.mesh = None
        self._synced = False
        super().__init__(src_ranks, trg_ranks)

    @property
    def synced(self):
        """bool: Tell if the DEC has already been synced."""
        return self._synced

    def synchronize(self):
        """Wrapper on DEC function."""
        if self._synced:
            return
        self._synced = True
        return super().synchronize()

    def attachLocalField(self, pfield):
        """Attach a local field"""
        self.pfield = pfield

    def sendData(self):
        """Send the field attached"""
        if self.isInSourceSide():
            self.sendToTarget(self.pfield.getField())
        else:
            self.sendToSource(self.pfield.getField())

    def recvData(self):
        """Received the field inside the field attached"""
        if self.isInSourceSide():
            field = self.receiveFromTarget()
        else:
            field = self.receiveFromSource()

        self.pfield.setArray(field.getArray())


class MEDCoupler:
    """Class handling the MEDCoupling related calls."""

    dec = log = None
    mesh_interf = mc_interf = mesh = None
    exch_fields = None
    debug = False

    __setattr__ = no_new_attributes(object.__setattr__)

    def supportOverlap(self):
        # To remove with IKDEC
        return hasattr(PMM, "InterpKernelDECWithOverlap")

    def __init__(self, logfunc=None, debug=False):
        self.dec = {}
        self.mesh_interf = self.mc_interf = self.mesh = None
        self.exch_fields = {}
        self.debug = debug
        self.log = logfunc if logfunc else logger

    def init_coupling(self, ranks1, ranks2, list_dec=None):
        """Start ParaMEDMEM coupling and DEC.

        Arguments:
            ranks1 (list[int]): List of ranks allocated to the first application.
            ranks2 (list[int]): List of ranks allocated to the second application.
            list_dec [list[MEDC.TypeOfField]] : list of dec to create
        """
        self.log("initializing ParaMEDMEM DEC")

        dec_to_create = list_dec
        if dec_to_create is None:
            dec_to_create = [MEDC.ON_CELLS, MEDC.ON_NODES, MEDC.ON_NODES_FE]

        if MEDC.ON_CELLS in dec_to_create:
            self.log(
                f"creating InterpKernelDEC on cells with " f"ranks1={ranks1} and ranks2={ranks2}",
                verbosity=2,
            )

            self.dec[MEDC.ON_CELLS] = ExtendedInterpKernelDECWithOverlap(ranks1, ranks2)
            self.dec[MEDC.ON_CELLS].setMethod("P0")

        if MEDC.ON_NODES in dec_to_create:
            self.log(
                f"creating InterpKernelDEC on nodes with " f"ranks1={ranks1} and ranks2={ranks2}",
                verbosity=2,
            )

            self.dec[MEDC.ON_NODES] = ExtendedInterpKernelDECWithOverlap(ranks1, ranks2)
            self.dec[MEDC.ON_NODES].setMethod("P1")

        if MEDC.ON_NODES_FE in dec_to_create:

            self.log(
                f"creating CFEMDEC on nodes with " f"ranks1={ranks1} and ranks2={ranks2}",
                verbosity=2,
            )

            self.dec[MEDC.ON_NODES_FE] = ExtendedCFEMDEC(ranks1, ranks2)

    def _create_paramesh(self, nodesIdsGlob, cellsIdsGlob):
        """Create the ParaMEDMEM mesh, support of coupling."""

        self.log("creating coupling mesh in memory", verbosity=2)
        for sup, dec in self.dec.items():
            if dec.isInSourceSide():
                group = dec.getSourceGrp()
            else:
                group = dec.getTargetGrp()
            dec.mesh = PMM.ParaMESH(self.mc_interf, group, "couplingMesh")
            dec.mesh.setCellGlobal(cellsIdsGlob)
            if sup != MEDC.ON_NODES:
                dec.mesh.setNodeGlobal(nodesIdsGlob)

    def _chech_mesh(self):
        assert self.mesh_interf.getNumberOfNodes() == self.mc_interf.getNumberOfNodes()
        assert self.mesh_interf.getNumberOfCells() == self.mc_interf.getNumberOfCells()

        connAster = self.mesh_interf.getConnectivity()

        for cellId in range(self.mesh_interf.getNumberOfCells()):
            nodesAst = connAster[cellId]
            nodesMec = self.mc_interf.getNodeIdsOfCell(cellId)

            if sorted(nodesAst) != sorted(nodesMec):
                raise RuntimeError(
                    f"Incompatible mesh connectivity of cell {cellId}: {nodesAst} vs {nodesMec}"
                )

        # TODO: add indirection and use it for medc <-> aster field conversion

    def create_mesh_interface(self, mesh, groupsOfCells):
        """Create the Medcoupling mesh of the interface.
           The mesh is restricted to a given list of groups of cells.

        Arguments:
            mesh (Mesh|ParallelMesh): mesh.
            groupsOfCells (list[str]): list of groups of cells.
        """

        self.log("creating interface mesh", verbosity=2)

        if MPI.ASTER_COMM_WORLD.Get_size() > 1:
            assert isinstance(mesh, ParallelMesh)

        self.mesh = mesh
        self.mesh_interf = mesh.restrict(groupsOfCells)

        groupsOfCells_res = []
        for grp in groupsOfCells:
            if self.mesh.hasGroupOfCells(grp, local=True):
                groupsOfCells_res.append(grp)

        mm = self.mesh_interf.createMedCouplingMesh()
        levels = mm.getGrpsNonEmptyLevels(groupsOfCells_res)
        assert len(levels) == 1, "Groups are not at one level"
        meshDimRelToMaxExt = levels[0]
        self.mc_interf = mm.getMeshAtLevel(meshDimRelToMaxExt)
        # check that each mesh has the same cell numbering
        self._chech_mesh()

        # add mesh and global ids - specific to CFEMDEC
        if self.mesh_interf.isParallel():
            nodesIdsGl = MEDC.DataArrayInt64(self.mesh_interf.getLocalToGlobalNodeIds())
            assert nodesIdsGl.getNbOfElems() == self.mesh_interf.getNumberOfNodes()
            cellsIdsGl = MEDC.DataArrayInt64(self.mesh_interf.getLocalToGlobalCellIds())
            assert cellsIdsGl.getNbOfElems() == self.mesh_interf.getNumberOfCells()

            # cell numbering could be incomplete
            ids = cellsIdsGl.findIdsStrictlyNegative()
            if ids.getNbOfElems() > 0:
                raise RuntimeError(
                    "Some cells have a negative Ids. This is forbidden. Update the global cell numbering"
                )
        else:
            nodesIdsGl = MEDC.DataArrayInt64(
                [i for i in range(self.mesh_interf.getNumberOfNodes())]
            )
            cellsIdsGl = MEDC.DataArrayInt64(
                [i for i in range(self.mesh_interf.getNumberOfCells())]
            )

        self._create_paramesh(nodesIdsGl, cellsIdsGl)

        if MEDC.ON_NODES_FE in self.dec:
            self.dec[MEDC.ON_NODES_FE].attachLocalMesh(self.mc_interf, nodesIdsGl)

    def get_field(self, name, silent=False):
        """Return a coupled field by name.

        Arguments:
            name (str): Field name.

        Returns:
            *Field*: pfield field or *None* if not found.
        """
        found = self.exch_fields.get(name)
        if found:
            return found
        if not silent and not found:
            msg = f"Field {name} was not defined beforehand!"
            self.log(msg)
            raise KeyError(msg)
        return found

    def restrict_field(self, field, cmps=[]):
        """Create a new field restricted to the interface mesh.

        Arguments:
            field (Field) aster field to restrict.

        Returns:
            SimpleField: restricted field.
        """

        loc = field.getLocalization()

        if loc == "NOEU":
            sfield = field.toSimpleFieldOnNodes()
        else:
            assert loc == "ELEM"
            sfield = field.toSimpleFieldOnCells()

        return sfield.transfert(self.mesh_interf, cmps)

    def extent_field(self, field, cmps=[]):
        """Create a new field extent to the whole mesh.

        Arguments:
            field (SimpleField) aster field to extent.

        Returns:
            SimpleField: extented field.
        """

        return field.transfert(self.mesh, cmps)

    def add_field(self, field_name, components, field_type):
        """Add a coupled field.

        Arguments:
            field_name (str): Field name.
            components (list[str]): Components of the field.
            field_type (str): On "NODES" or "NODES_FE" or "CELLS".
        """
        if not self.get_field(field_name, silent=True):
            assert field_type in ("NODES", "NODES_FE", "CELLS")
            conv = {"NODES": MEDC.ON_NODES, "NODES_FE": MEDC.ON_NODES_FE, "CELLS": MEDC.ON_CELLS}
            sup = conv[field_type]
            dec = self.dec[sup]
            if field_type == "CELLS":
                nature = MEDC.IntensiveConservation
            else:
                nature = MEDC.IntensiveMaximum

            topo = PMM.ComponentTopology(len(components))
            pfield = CoupledField(sup, MEDC.ONE_TIME, dec, topo)
            field = pfield.getField()
            field.setName(field_name)
            pfield.setNature(nature)
            pfield.fillWithZero()
            pfield.setInfoOnComponents(components)

            self.exch_fields[field_name] = pfield

            self.log(f"add field {field_name!r} on {field_type}...")
            self.log(repr(pfield.getArray()), verbosity=2)

    def send(self, fields):
        """Send fields to the partner code with ParaMEDMEM.

        Arguments:
            fields (dict[*ParaFIELD*]): Fields to send.
        """
        for field_name, field in fields.items():
            pfield = self.get_field(field_name)
            dec = pfield.dec
            sup_name = pfield.sup_name
            self.log(f"sending field {field_name!r} on {sup_name}...")
            # update values
            pfield.setArray(field.getArray())
            if self.debug:
                self.log(repr(field), verbosity=2)
            dec.attachLocalField(pfield)
            self.log("sync...", verbosity=2)
            dec.synchronize()
            start = time.perf_counter()
            self.log("sendData...", verbosity=2)
            dec.sendData()
            delta = time.perf_counter() - start
            self.log(f"in {delta} sec.")

        self.log("pmm_send: done", verbosity=2)

    def recv(self, fields_names):
        """Receive fields from the partner code with ParaMEDMEM.

        Arguments:
            fields_names (list[str]): Fields names.

        Returns:
            dict[*CoupledField*]: Received fields.
        """
        fields = {}

        for field_name in fields_names:
            pfield = self.get_field(field_name)
            dec = pfield.dec
            sup_name = pfield.sup_name
            self.log(f"waiting for field {field_name!r} on {sup_name}...")
            dec.attachLocalField(pfield)
            fields[field_name] = pfield
            self.log("sync...", verbosity=2)
            dec.synchronize()
            self.log("recvData...", verbosity=2)
            start = time.perf_counter()
            dec.recvData()
            delta = time.perf_counter() - start
            self.log(f"received field {field_name!r} on {sup_name} in {delta} sec...", verbosity=2)
            if self.debug:
                self.log(repr(pfield.getField()), verbosity=2)

        self.log("pmm_recv: done", verbosity=2)
        return fields

    def _medcfield2aster(self, mc_field):
        """Convert MEDCouplingField to FieldOnNodes/Cells

        Arguments:
            mc_field (*MEDCouplingField*): MEDCoupling field.

        Returns:
            *SimpleField*: aster field.
        """
        if mc_field.getTypeOfField() in (MEDC.ON_NODES, MEDC.ON_NODES_FE):
            if isinstance(mc_field, CoupledField):
                sfield = SimpleFieldOnNodesReal.fromMedCouplingField(
                    mc_field.getField(), self.mesh_interf
                )
            else:
                sfield = SimpleFieldOnNodesReal.fromMedCouplingField(mc_field, self.mesh_interf)
        else:
            if isinstance(mc_field, CoupledField):
                sfield = SimpleFieldOnCellsReal.fromMedCouplingField(
                    mc_field.getField(), self.mesh_interf
                )
            else:
                sfield = SimpleFieldOnCellsReal.fromMedCouplingField(mc_field, self.mesh_interf)

        return sfield

    def import_field(self, mc_field, physq, symbname, model=None):
        """Convert a MEDCoupling field defined on the interface as
        a code_aster field defined on the whole mesh.

        Arguments:
            mc_field (*MEDCouplingField*): MEDCoupling field.
            physq (str): Physical quantity of field (e.g. DEPL_R).
            symbname (str): Symbolic name of field (e.g. SIEQ_ELGA).
            model (Model): model to use for FieldOnCells.

        Returns:
            *FieldOnNodesReal*: code_aster field defined on the whole mesh.
        """

        internal_desc = "-".join((physq, symbname))
        mc_field.setDescription(internal_desc)

        field = self._medcfield2aster(mc_field)

        if mc_field.getTypeOfField() in (MEDC.ON_NODES_FE, MEDC.ON_NODES):
            return self.extent_field(field).toFieldOnNodes()
        else:
            fed = model.getFiniteElementDescriptor().restrict(self.mesh_interf.getGroupsOfCells())
            return self.extent_field(field).toFieldOnCells(fed)

    def export_field(self, field, field_name="COUPLINGFIELD", cmps=[]):
        """Convert a code_aster field defined on the whole mesh to
            a MEDCoupling field defined on the interface.

        Arguments:
            field *FieldOnNodesReal*: code_aster field defined on the whole mesh.
            field_name (str): name of the field (like `DEPL`) (default: field's name)
            cmps (list[str]): list of components. (default: all)

        Returns:
            *MEDCouplingField*: MEDCoupling field.
        """

        if len(cmps) == 0:
            cmps = field.getComponents()

        assert field.getMesh() == self.mesh

        field_interf = self.restrict_field(field, cmps)
        pfield = field_interf.toMedCouplingField(self.mc_interf, field_name)

        return pfield

    def import_displacement(self, mc_displ):
        """Convert a MEDCoupling displacement field defined on the interface as
        a code_aster field.

        Arguments:
            mc_displ (*MEDCouplingField*): MEDCoupling displacement field.

        Returns:
            FieldOnNodesReal: code_aster displacement field.
        """

        return self.import_field(mc_displ, "DEPL_R", "DEPL")

    def export_displacement(self, displ, field_name="DEPL"):
        """Create a MEDCoupling field of displacement reduced on the interface mesh.

        Arguments:
            displ (FieldOnNodesReal): code_aster displacement field.
            field_name (str): Field name. (default: field's name)

        Returns:
            *MEDCouplingFieldDouble*: Displacement field.
        """

        return self.export_field(displ, field_name, ["DX", "DY", "DZ"])

    def import_velocity(self, mc_velo):
        """Convert a MEDCoupling velocity field defined on the interface as
        a code_aster field.

        Arguments:
            mc_velo (*MEDCouplingField*): MEDCoupling velocity field.

        Returns:
            FieldOnNodesReal: code_aster velocity field.
        """

        return self.import_displacement(mc_velo, "DEPL_R", "VITE")

    def export_velocity(self, velo, field_name="VELOCITY"):
        """Create a MEDCoupling field of velocity reduced on the interface mesh.

        Arguments:
            velo (FieldOnNodesReal): code_aster velocity field.
            field_name (str): Field name. (default: field's name)

        Returns:
            *MEDCouplingFieldDouble*: Velocity field.
        """

        return self.export_displacement(velo, field_name)

    def export_temperature(self, temp, field_name="TEMP"):
        """Create a MEDCoupling field of temperature reduced on the interface mesh.

        Arguments:
            temp (FieldOnNodesReal): code_aster thermal field.
            field_name (str): Field name. (default: field's name)

        Returns:
            *MEDCouplingFieldDouble*: Thermal field on cells.
        """

        return self.export_field(temp, field_name, ["TEMP"])

    def import_temperature(self, mc_temp):
        """Convert a MEDCoupling thermal field as a code_aster field.

        Arguments:
            mc_temp (*MEDCouplingFieldDouble*): MEDCoupling thermal field.

        Returns:
            FieldOnNodesReal: code_aster thermal field.
        """

        return self.import_field(mc_temp, "TEMP_R", "TEMP")

    def export_pressure(self, pres, field_name="PRES"):
        """Create a MEDCoupling field of pressure reduced on the interface mesh.

        Arguments:
            press (*FieldOnNodesReal*): code_aster pressure field.
            field_name (str): Field name. (default: field's name)

        Returns:
            *MEDCouplingFieldDouble*: Pressure field on cells.
        """

        return self.export_field(pres, field_name, ["PRES"])

    def import_pressure(self, mc_pres):
        """Convert a MEDCoupling pressure field as a code_aster field.

        Arguments:
            mc_pres (*MEDCouplingFieldDouble*): MEDCoupling pressure field.

        Returns:
            *FieldOnNodesReal*: code_aster pressure field.
        """

        return self.import_field(mc_pres, "PRES_R", "PRES")

    def import_fluidforces(self, mc_fluidf, model, time=0.0):
        """Convert a MEDCoupling fluid forces field as a code_aster field.

        Arguments:
            mc_fluidf (*MEDCouplingField*): MEDCoupling fluid forces field.
            model (Model): Mechanical model.
            time (float): Time of assignment.

        Returns:
            *LoadResult*: surface forces load.
        """

        forc_elem = self.import_field(mc_fluidf, "FORC_R", "FORC", model)

        evol_char = LoadResult()
        evol_char.allocate(1)
        evol_char.setModel(model, 0)
        evol_char.setTime(time, 0)
        evol_char.setField(forc_elem, "FSUR_3D", 0)
        evol_char.build()

        return evol_char
