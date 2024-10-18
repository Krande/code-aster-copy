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

"""
Definition of a convenient object to synchronize MEDCoupling fields.
"""

import os

import medcoupling as MEDC
import ParaMEDMEM as PMM

from ..Commands import CREA_RESU, LIRE_CHAMP, PROJ_CHAMP
from ..Objects import Mesh
from ..Utilities.logger import logger
from ..Utilities.mpi_utils import MPI


class CoupledField:
    """Define the properties of an exchanged field.

    Attributes:
        name (str): Field name.
        components (list[str]): Names of the components.
        support (int): Support: ON_NODES or ON_CELLS.
        size (int): Number of values to be exchanged on the interface.
        para_field (*ParaFIELD*): *ParaMEDMEM.ParaFIELD* object.
    """

    def __init__(self, name, components, support, size, para_field):
        self.name = name
        self.components = components
        assert support in (MEDC.ON_NODES, MEDC.ON_CELLS), support
        self.support = support
        self.size = size
        self.para_field = para_field


class ExtendedDEC:
    """Object that represents a DEC and the necessary properties.

    Arguments:
        mcdec (*PMM.InterpKernelDEC*): ParaMEDMEM InterpKernelDEC object.
    """

    def __init__(self, mcdec):
        self.dec = mcdec
        self.mesh = None
        self.interf_size = 0
        self._synced = False

    def isInSourceSide(self):
        """Wrapper on DEC function."""
        return self.dec.isInSourceSide()

    def getSourceGrp(self):
        """Wrapper on DEC function."""
        return self.dec.getSourceGrp()

    def getTargetGrp(self):
        """Wrapper on DEC function."""
        return self.dec.getTargetGrp()

    def attachLocalField(self, field):
        """Wrapper on DEC function."""
        return self.dec.attachLocalField(field)

    def setMethod(self, method):
        """Wrapper on DEC function."""
        return self.dec.setMethod(method)

    @property
    def synced(self):
        """bool: Tell if the DEC has already been synced."""
        return self._synced

    def synchronize(self):
        """Wrapper on DEC function."""
        if self._synced:
            return
        self._synced = True
        return self.dec.synchronize()

    def sendData(self):
        """Wrapper on DEC function."""
        return self.dec.sendData()

    def recvData(self):
        """Wrapper on DEC function."""
        return self.dec.recvData()


class MEDCoupler:
    """Class handling the MEDCoupling related calls."""

    def __init__(self, logfunc=None):
        self.dec = {MEDC.ON_NODES: {}, MEDC.ON_CELLS: {}}
        self.xdec = {MEDC.ON_NODES: {}, MEDC.ON_CELLS: {}}
        self.interf_mesh = self.interf_mc = self.mesh = None
        self.matr_proj = None
        self.exch_fields = {}

        self.log = logfunc if logfunc else logger

    def init_coupling(self, fields_name, ranks1, ranks2):
        """Start ParaMEDMEM coupling and DEC.

        Arguments:
            fields_name (list[str]): list of field names.
            ranks1 (list[int]): List of ranks allocated to the first application.
            ranks2 (list[int]): List of ranks allocated to the second application.
        """
        self.log("initializing ParaMEDMEM InterpKernelDEC")
        # Creating the send and recv DEC's
        for sup in (MEDC.ON_CELLS, MEDC.ON_NODES):
            sup_name = "nodes" if sup == MEDC.ON_NODES else "cells"
            key = f"{ranks1},{ranks2}"
            if key not in self.dec[sup]:
                self.log(
                    f"creating InterpKernelDEC on {sup_name} with "
                    f"ranks1={ranks1} and ranks2={ranks2}",
                    verbosity=2,
                )
                self.dec[sup][key] = PMM.InterpKernelDEC(ranks1, ranks2)
            for name in fields_name:
                self.xdec[sup][name] = ExtendedDEC(self.dec[sup][key])
                self.xdec[sup][name].setMethod("P0" if sup == MEDC.ON_CELLS else "P1")

    def _create_paramesh(self):
        """Create the ParaMEDMEM mesh, support of coupling."""

        self.log("creating coupling mesh in memory", verbosity=2)
        for sup in (MEDC.ON_CELLS, MEDC.ON_NODES):
            for name in self.xdec[sup]:
                dec = self.xdec[sup][name]
                if dec.isInSourceSide():
                    group = dec.getSourceGrp()
                else:
                    group = dec.getTargetGrp()
                dec.mesh = PMM.ParaMESH(self.interf_mc, group, "couplingMesh")

    def create_mesh_interface(self, mesh, groupsOfCells):
        """Create the Medcoupling mesh of the interface.
           The mesh is restricted to a given list of groups of cells.

        Arguments:
            mesh (Mesh|ParallelMesh): mesh.
            groupsOfCells (list[str]): list of groups of cells.
        """

        self.log("creating interface mesh", verbosity=2)

        self.mesh = mesh

        mm = self.mesh.createMedCouplingMesh()
        levels = mm.getGrpsNonEmptyLevels(groupsOfCells)
        assert len(levels) == 1, "Groups are not at one level"
        self.meshDimRelToMaxExt = levels[0]
        self.interf_ids = mm.getGroupsArr(self.meshDimRelToMaxExt, groupsOfCells)
        self.interf_ids.setName("interf")
        self.interf_mc = mm.getMeshAtLevel(0)[self.interf_ids]
        self.interf_mc.setName("interface")

        meshMEDFile = MEDC.MEDFileUMesh()
        meshMEDFile.setMeshAtLevel(0, self.interf_mc)
        self.mesh_interf = Mesh()
        self.mesh_interf.buildFromMedCouplingMesh(meshMEDFile)

        self._create_paramesh()

    def get_field(self, name, silent=False):
        """Return a coupled field by name.

        Arguments:
            name (str): Field name.

        Returns:
            *Field*: Exchanged field or *None* if not found.
        """
        found = self.exch_fields.get(name)
        if found:
            return found
        if not silent and not found:
            msg = f"Field {name} was not defined beforehand!"
            self.log(msg)
            raise KeyError(msg)
        return found

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
            dec = self.xdec[sup][field_name]
            if field_type == "CELLS":
                nature = MEDC.IntensiveConservation
                nb_tuples = self.interf_mc.getNumberOfCells()
            else:
                nature = MEDC.IntensiveMaximum
                nb_tuples = self.interf_mc.getNumberOfNodes()

            topo = PMM.ComponentTopology(len(components))
            pfield = PMM.ParaFIELD(sup, MEDC.ONE_TIME, dec.mesh, topo)
            pfield.getField().setName(field_name)
            pfield.getField().setNature(nature)
            pfield.getField().getArray().setInfoOnComponents(components)

            if not dec.interf_size:
                dec.interf_size = nb_tuples * len(components)
            assert (
                nb_tuples * len(components) == dec.interf_size
            ), "inconsistent number of components: %i vs. %i" % (
                nb_tuples * len(components),
                dec.interf_size,
            )

            field = CoupledField(field_name, components, sup, dec.interf_size, pfield)
            self.exch_fields[field_name] = field

    def sync(self):
        """Synchronize the parallel DEC objects."""
        for sup in (MEDC.ON_CELLS, MEDC.ON_NODES):
            for name in self.xdec[sup]:
                dec = self.xdec[sup][name]
                if dec.synced:
                    continue
                cfield = [field for field in self.exch_fields.values() if field.support == sup]
                if len(cfield) > 0:
                    cfield = cfield[0]
                    pfield = cfield.para_field.getField()
                    dec.attachLocalField(pfield)
                    sup_name = "nodes" if sup == MEDC.ON_NODES else "cells"
                    self.log(
                        f"synchronize DEC #{sup}-{name} for field {name} on {sup_name}, "
                        f"interface size: {dec.interf_size}",
                        verbosity=2,
                    )
                    dec.synchronize()

    def send(self, fields):
        """Send fields to the partner code with ParaMEDMEM.

        Arguments:
            fields (dict[*ParaFIELD*]): Fields to send.
        """
        if not fields:
            return
        support = None
        for field_name in fields:
            field = fields[field_name]
            exchanged = self.get_field(field_name)
            if support is None:
                support = exchanged.support
            assert support == support, "all fields must be on the same support type."
            sup_name = "nodes" if support == MEDC.ON_NODES else "cells"
            self.log(f"sending field {field_name!r} on {sup_name}...")
            pfield = exchanged.para_field.getField()
            array = field.getArray()
            pfield.setArray(array)
            pfield.setName(field_name)
            self.log(repr(pfield), verbosity=2)
            # self.log(array, verbosity=2)
            dec = self.xdec[support][field_name]
            dec.attachLocalField(pfield)
            dec.synchronize()
            self.log(f"sendData...")
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
        if not fields_names:
            return fields
        support = None
        for field_name in fields_names:
            exchanged = self.get_field(field_name)
            if support is None:
                support = exchanged.support
            assert support == support, "all fields must be on the same support type."
            sup_name = "nodes" if support == MEDC.ON_NODES else "cells"
            self.log(f"waiting for field {field_name!r} on {sup_name}...")
            pfield = exchanged.para_field.getField()
            pfield.setValues([0.0] * exchanged.size)
            array = pfield.getArray()
            array.rearrange(len(exchanged.components))
            array.setInfoOnComponents(exchanged.components)
            dec = self.xdec[support][field_name]
            dec.attachLocalField(pfield)
            fields[field_name] = pfield
            dec.synchronize()
            self.log(f"recvData...", verbosity=2)
            dec.recvData()
            self.log(repr(pfield), verbosity=2)
            # self.log(array, verbosity=2)
        self.log(f"pmm_recv: done", verbosity=2)
        return fields

    @staticmethod
    def _write_field2med(field, filename):
        """Write field on disk using MED format.

        Arguments:
            field (FieldOn*): aster field.
            filename (str): name of MED file.
        """
        if os.path.exists(filename):
            os.remove(filename)
        logger.debug("writing file {0!r}...".format(filename), flush=True)
        field.printMedFile(filename)

    def _medcfield2aster(self, mc_field, field_type, time=0.0):
        """Convert MEDCouplingField to FieldOnNodes/Cells

        Arguments:
            mc_field (*MEDCouplingField*): MEDCoupling field.
            field_type (str): type of the field (like `NOEU_DEPL_R`)
            time (float, optional): Time of assignment.

        Returns:
            *Field*: aster field.
        """

        mc_mesh = mc_field.getMesh()
        tmpfile = "fort.77"
        MEDC.WriteUMesh(tmpfile, mc_mesh, True)
        MEDC.WriteFieldUsingAlreadyWrittenMesh(tmpfile, mc_field)
        field = LIRE_CHAMP(
            UNITE=77,
            MAILLAGE=self.mesh_interf,
            PROL_ZERO="OUI",
            NOM_MED=mc_field.getName(),
            TYPE_CHAM=field_type,
            NOM_CMP_IDEM="OUI",
            INST=time,
            INFO=1,
        )
        os.remove(tmpfile)
        return field

    def import_field(self, mc_field, field_type, time=0.0):
        """Convert a MEDCoupling field defined on the interface as
        a code_aster field defined on the whole mesh.

        Arguments:
            mc_field (*MEDCouplingField*): MEDCoupling field.
            field_type (str): type of the field (like `NOEU_DEPL_R`)
            time (float, optional): Time of assignment.

        Returns:
            *FieldOnNodes*: code_aster field defined on the whole mesh.
        """

        field = self._medcfield2aster(mc_field, field_type, time)
        return self.project_field(field)

    def import_displacement(self, mc_displ, time=0.0):
        """Convert a MEDCoupling displacement field defined on the interface as
        a code_aster field.

        Arguments:
            mc_displ (*MEDCouplingField*): MEDCoupling displacement field.
            time (float, optional): Time of assignment.

        Returns:
            *FieldOnNodesReal*: code_aster displacement field.
        """

        return self.import_field(mc_displ, "NOEU_DEPL_R", time)

    def export_displacement(self, displ, field_name):
        """Create a MEDCoupling field of displacement reduced on the interface mesh.

        Arguments:
            displ (*FieldOnNodes*): code_aster displacement field.
            field_name (str): Field name.

        Returns:
            *MEDCouplingField*: Displacement field.
        """
        cmps = [cmp for cmp in displ.getComponents() if cmp in ["DX", "DY", "DZ"]]
        if len(cmps) != len(displ.getComponents()):
            displ = displ.restrict(cmps=cmps)

        filename = "displ_%i.med" % (MPI.ASTER_COMM_WORLD.rank)
        self._write_field2med(displ, filename)

        fieldname = displ.getName()[:8]
        mc_displ = MEDC.ReadFieldNode(
            filename, displ.getMesh().getName(), self.meshDimRelToMaxExt, fieldname, -1, -1
        )

        displ = mc_displ[self.interf_ids]
        displ.setName(field_name)
        displ.setNature(MEDC.IntensiveMaximum)
        displ.setMesh(self.interf_mc)
        array = displ.getArray()
        array.setInfoOnComponents(cmps)
        displ.checkConsistencyLight()
        os.remove(filename)

        return displ

    def export_temperature(self, temp, field_name):
        """Create a MEDCoupling field of temperature reduced on the interface mesh.

        Arguments:
            temp (*FieldOnNodes*): code_aster thermal field.
            field_name (str): Field name.

        Returns:
            *MEDCouplingField*: Pressure field on cells.
        """

        filename = "temp_%i.med" % (MPI.ASTER_COMM_WORLD.rank)
        self._write_field2med(temp, filename)

        fieldname = temp.getName()[:8]
        mc_temp = MEDC.ReadFieldNode(
            filename, temp.getMesh().getName(), self.meshDimRelToMaxExt, fieldname, -1, -1
        )

        temp = mc_temp[self.interf_ids]
        temp.setName(field_name)
        temp.setNature(MEDC.IntensiveMaximum)
        temp.setMesh(self.interf_mc)
        temp.getArray().setInfoOnComponents(["TEMP"])
        temp.checkConsistencyLight()
        os.remove(filename)

        return temp

    def import_temperature(self, mc_temp, time=0.0):
        """Convert a MEDCoupling pemperature field as a code_aster field.

        Arguments:
            mc_temp (*MEDCouplingField*): MEDCoupling temp field.
            time (float, optional): Time of assignment.

        Returns:
            *FieldOnNodes*: code_aster thermal field.
        """

        return self.import_field(mc_temp, "NOEU_TEMP_R", time)

    def import_fluidforces(self, mc_fluidf, field_name, time=0.0):
        """Convert a MEDCoupling pressure field as a code_aster field.

        Arguments:
            mc_fluidf (*MEDCouplingField*): MEDCoupling pressure field.
            field_name (str): Field name.
            time (float, optional): Time of assignment.

        Returns:
            *LoadResult*: code_aster pressure as *LoadResult*.
        """

        forces = self.import_field(mc_fluidf, "ELEM_FORC_R", time)

        forces = CREA_RESU(
            TYPE_RESU="EVOL_CHAR",
            OPERATION="AFFE",
            AFFE=_F(NOM_CHAM="FORC_NODA", CHAM_GD=forces, MODELE=self.model_interf, INST=time),
        )
        # projection on the entire model
        proj_forces = PROJ_CHAMP(
            METHODE="COLLOCATION",
            RESULTAT=forces,
            MODELE_1=self.model_interf,
            MODELE_2=self.model,
            # VIS_A_VIS=_F(
            #     GROUP_MA_1="interf",
            #     GROUP_MA_2="interf",
            # ),
            TOUT_ORDRE="OUI",
            PROL_ZERO="OUI",
        )
        return proj_forces

    def project_field(self, field):
        """Project field from ther interface to the whole mesh and
        extend by zero where the field does not exist.

        Arguments:
            field (FieldOn***): field defined on the interface mesh.

        Returns
            (FieldOn***): field projected on the whole mesh.
        """

        if self.matr_proj is None:
            self.matr_proj = PROJ_CHAMP(
                METHODE="COLLOCATION",
                PROJECTION="NON",
                MAILLAGE_1=self.mesh_interf,
                MAILLAGE_2=self.mesh,
            )

        return PROJ_CHAMP(CHAM_GD=field, MATR_PROJECTION=self.matr_proj, PROL_ZERO="OUI")
