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

"""
Definition of a convenient object to synchronize MEDCoupling fields.
"""

from ..Objects import LoadResult, SimpleFieldOnCellsReal, SimpleFieldOnNodesReal
from ..Utilities import ParaMEDMEM as PMM
from ..Utilities import logger, no_new_attributes
from ..Utilities import medcoupling as MEDC


class CoupledField(PMM.ParaFIELD):
    """Define the properties of an coupled field.

    Attributes:
        support (TypeOfField): Support: ON_NODES or ON_CELLS.
        td (TypeOfTimeDiscretization): time discretization
        dec (ExtendedDEC): dec.
        topo (ComponentTopology): topology.
    """

    def __init__(self, sup, td, dec, topo):
        assert sup in (MEDC.ON_NODES, MEDC.ON_CELLS), sup
        self.dec = dec
        super().__init__(sup, td, dec.mesh, topo)

    def fillWithZero(self):
        self.getField().getArray().fillWithZero()

    def setArray(self, array):
        assert self.getNbOfElems() == array.getNbOfElems()
        self.getField().setArray(array)

    def getNbOfElems(self):
        return self.getField().getArray().getNbOfElems()

    @property
    def sup(self):
        return self.getField().getTypeOfField()


class ExtendedDEC(PMM.InterpKernelDEC):
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


class MEDCoupler:
    """Class handling the MEDCoupling related calls."""

    dec = log = None
    mesh_interf = mc_interf = mesh = None
    exch_fields = None
    debug = False

    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, logfunc=None, debug=False):
        self.dec = {MEDC.ON_NODES: None, MEDC.ON_CELLS: None}
        self.mesh_interf = self.mc_interf = self.mesh = None
        self.exch_fields = {}
        self.debug = debug
        self.log = logfunc if logfunc else logger

    def init_coupling(self, ranks1, ranks2):
        """Start ParaMEDMEM coupling and DEC.

        Arguments:
            ranks1 (list[int]): List of ranks allocated to the first application.
            ranks2 (list[int]): List of ranks allocated to the second application.
        """
        self.log("initializing ParaMEDMEM InterpKernelDEC")
        # Creating the send and recv DEC's
        for sup in (MEDC.ON_CELLS, MEDC.ON_NODES):
            sup_name = "nodes" if sup == MEDC.ON_NODES else "cells"
            self.log(
                f"creating InterpKernelDEC on {sup_name} with "
                f"ranks1={ranks1} and ranks2={ranks2}",
                verbosity=2,
            )
            self.dec[sup] = ExtendedDEC(ranks1, ranks2)
            self.dec[sup].setMethod("P0" if sup == MEDC.ON_CELLS else "P1")

    def _create_paramesh(self):
        """Create the ParaMEDMEM mesh, support of coupling."""

        self.log("creating coupling mesh in memory", verbosity=2)
        for dec in self.dec.values():
            if dec.isInSourceSide():
                group = dec.getSourceGrp()
            else:
                group = dec.getTargetGrp()
            dec.mesh = PMM.ParaMESH(self.mc_interf, group, "couplingMesh")

    def create_mesh_interface(self, mesh, groupsOfCells):
        """Create the Medcoupling mesh of the interface.
           The mesh is restricted to a given list of groups of cells.

        Arguments:
            mesh (Mesh|ParallelMesh): mesh.
            groupsOfCells (list[str]): list of groups of cells.
        """

        self.log("creating interface mesh", verbosity=2)

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

        self._create_paramesh()

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
            field_type (str): On "NODES" or "CELLS".
        """
        if not self.get_field(field_name, silent=True):
            assert field_type in ("NODES", "CELLS")
            sup = MEDC.ON_CELLS if field_type == "CELLS" else MEDC.ON_NODES
            dec = self.dec[sup]
            if field_type == "CELLS":
                nature = MEDC.IntensiveConservation
            else:
                nature = MEDC.IntensiveMaximum

            topo = PMM.ComponentTopology(len(components))
            pfield = CoupledField(sup, MEDC.ONE_TIME, dec, topo)
            field = pfield.getField()
            field.setName(field_name)
            field.setNature(nature)
            array = field.getArray()
            array.fillWithZero()
            array.setInfoOnComponents(components)

            self.exch_fields[field_name] = pfield

            self.log(f"add field {field_name!r} on {field_type}...")
            self.log(repr(array), verbosity=2)

    def send(self, fields):
        """Send fields to the partner code with ParaMEDMEM.

        Arguments:
            fields (dict[*ParaFIELD*]): Fields to send.
        """
        for field_name, field in fields.items():
            pfield = self.get_field(field_name)
            dec = pfield.dec
            sup_name = "nodes" if pfield.sup == MEDC.ON_NODES else "cells"
            self.log(f"sending field {field_name!r} on {sup_name}...")
            # update values
            pfield.setArray(field.getArray())
            if self.debug:
                self.log(repr(field), verbosity=2)
            dec.attachLocalField(pfield)
            self.log("sync...", verbosity=2)
            dec.synchronize()
            self.log("sendData...")
            dec.sendData()
        self.log("pmm_send: done", verbosity=2)

    def recv(self, fields_names):
        """Receive fields from the partner code with ParaMEDMEM.

        Arguments:
            fields_names (list[str]): Fields names.

        Returns:
            dict[*ParaFIELD*]: Received fields.
        """
        fields = {}

        for field_name in fields_names:
            pfield = self.get_field(field_name)
            dec = pfield.dec
            sup_name = "nodes" if pfield.sup == MEDC.ON_NODES else "cells"
            self.log(f"waiting for field {field_name!r} on {sup_name}...")
            dec.attachLocalField(pfield)
            fields[field_name] = pfield.getField()
            self.log("sync...", verbosity=2)
            dec.synchronize()
            self.log("recvData...", verbosity=2)
            dec.recvData()
            self.log(f"received field {field_name!r} on {sup_name}...")
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
        if mc_field.getTypeOfField() == MEDC.ON_NODES:
            sfield = SimpleFieldOnNodesReal.fromMedCouplingField(mc_field, self.mesh_interf)
        else:
            sfield = SimpleFieldOnCellsReal.fromMedCouplingField(mc_field, self.mesh_interf)

        return sfield

    def import_field(self, mc_field, model=None):
        """Convert a MEDCoupling field defined on the interface as
        a code_aster field defined on the whole mesh.

        Arguments:
            mc_field (*MEDCouplingField*): MEDCoupling field.
            model (Model): model to use for FieldOnCells.

        Returns:
            *FieldOnNodesReal*: code_aster field defined on the whole mesh.
        """

        field = self._medcfield2aster(mc_field)

        if mc_field.getTypeOfField() == MEDC.ON_NODES:
            return self.extent_field(field).toFieldOnNodes()
        else:
            fed = model.getFiniteElementDescriptor().restrict(self.mesh_interf.getGroupsOfCells())
            return self.extent_field(field).toFieldOnCells(fed)

    def export_field(self, field, field_name=None, cmps=[]):
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
        pfield = field_interf.toMedCouplingField(self.mc_interf)
        if field_name:
            pfield.setName(field_name)

        return pfield

    def import_displacement(self, mc_displ):
        """Convert a MEDCoupling displacement field defined on the interface as
        a code_aster field.

        Arguments:
            mc_displ (*MEDCouplingField*): MEDCoupling displacement field.

        Returns:
            FieldOnNodesReal: code_aster displacement field.
        """

        mc_displ.getArray().setName("DEPL_R")

        return self.import_field(mc_displ)

    def export_displacement(self, displ, field_name=None):
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
        return self.import_displacement(mc_velo)

    def export_velocity(self, velo, field_name=None):
        """Create a MEDCoupling field of velocity reduced on the interface mesh.

        Arguments:
            velo (FieldOnNodesReal): code_aster velocity field.
            field_name (str): Field name. (default: field's name)

        Returns:
            *MEDCouplingFieldDouble*: Velocity field.
        """

        return self.export_displacement(velo, field_name)

    def export_temperature(self, temp, field_name=None):
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

        mc_temp.getArray().setName("TEMP_R")

        return self.import_field(mc_temp)

    def export_pressure(self, pres, field_name=None):
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

        mc_pres.getArray().setName("PRES_R")

        return self.import_field(mc_pres)

    def import_fluidforces(self, mc_fluidf, model, time=0.0):
        """Convert a MEDCoupling fluid forces field as a code_aster field.

        Arguments:
            mc_fluidf (*MEDCouplingField*): MEDCoupling fluid forces field.
            model (Model): Mechanical model.
            time (float): Time of assignment.

        Returns:
            *LoadResult*: surface forces load.
        """

        mc_fluidf.getArray().setName("FORC_R")
        forc_elem = self.import_field(mc_fluidf, model)

        evol_char = LoadResult()
        evol_char.allocate(1)
        evol_char.setModel(model, 0)
        evol_char.setTime(time, 0)
        evol_char.setField(forc_elem, "FSUR_3D", 0)
        evol_char.build()

        return evol_char
