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

from ..Commands import LIRE_CHAMP, PROJ_CHAMP, AFFE_MODELE
from ..Objects import Mesh, LoadResult
from ..Utilities import logger, no_new_attributes
from ..Utilities.mpi_utils import MPI


class CoupledField:
    """Define the properties of an exchanged field.

    Attributes:
        name (str): Field name.
        components (list[str]): Names of the components.
        support (int): Support: ON_NODES or ON_CELLS.
        dec (ExtendedDEC): dec.
        size (int): Number of values to be exchanged on the interface.
        para_field (*ParaFIELD*): *ParaMEDMEM.ParaFIELD* object.
    """

    def __init__(self, name, components, support, dec, size, para_field):
        self.name = name
        self.components = components
        assert support in (MEDC.ON_NODES, MEDC.ON_CELLS), support
        self.support = support
        self.size = size
        self.para_field = para_field
        self.dec = dec


class ExtendedDEC:
    """Object that represents a DEC and the necessary properties.

    Arguments:
        mcdec (*PMM.InterpKernelDEC*): ParaMEDMEM InterpKernelDEC object.
    """

    dec = mesh = interf_size = _synced = None
    __setattr__ = no_new_attributes(object.__setattr__)

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

    dec = None
    mesh_interf = interf_mc = mesh = interf_ids = None
    model_interf = None
    matr_proj = meshDimRelToMaxExt = None
    exch_fields = None
    log = None

    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, logfunc=None):
        self.dec = {MEDC.ON_NODES: None, MEDC.ON_CELLS: None}
        self.mesh_interf = self.interf_mc = self.mesh = None
        self.matr_proj = None
        self.exch_fields = {}

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
            self.dec[sup] = ExtendedDEC(PMM.InterpKernelDEC(ranks1, ranks2))
            self.dec[sup].setMethod("P0" if sup == MEDC.ON_CELLS else "P1")

    def _create_paramesh(self):
        """Create the ParaMEDMEM mesh, support of coupling."""

        self.log("creating coupling mesh in memory", verbosity=2)
        for dec in self.dec.values():
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

        groupsOfCells_res = []
        for grp in groupsOfCells:
            if self.mesh.hasGroupOfCells(grp, local=True):
                groupsOfCells_res.append(grp)

        mm = self.mesh.createMedCouplingMesh()
        levels = mm.getGrpsNonEmptyLevels(groupsOfCells_res)
        assert len(levels) == 1, "Groups are not at one level"
        self.meshDimRelToMaxExt = levels[0]
        self.interf_ids = mm.getGroupsArr(self.meshDimRelToMaxExt, groupsOfCells_res)
        self.interf_ids.setName("interf")
        self.interf_mc = mm.getMeshAtLevel(self.meshDimRelToMaxExt)[self.interf_ids]
        self.interf_mc.setName("interface")
        # remove orphelan nodes
        self.interf_mc.zipCoords()

        meshMEDFile = MEDC.MEDFileUMesh()
        meshMEDFile.setMeshAtLevel(self.meshDimRelToMaxExt, self.interf_mc)

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
            dec = self.dec[sup]
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

            field = CoupledField(
                field_name, components, sup, dec, nb_tuples * len(components), pfield
            )
            self.exch_fields[field_name] = field

    def sync(self):
        """Synchronize the parallel DEC objects."""
        for sup, dec in self.dec.items():
            if dec.synced:
                continue
            cfield = [field for field in self.exch_fields.values() if field.support == sup]
            if len(cfield) > 0:
                cfield = cfield[0]
                pfield = cfield.para_field.getField()
                dec.attachLocalField(pfield)
                sup_name = "nodes" if sup == MEDC.ON_NODES else "cells"
                self.log(
                    f"synchronize DEC #{sup}  on {sup_name}, ",
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
        for sup, dec in self.dec.items():
            nb_field = 0
            for field_name, field in fields.items():
                exchanged = self.get_field(field_name)
                support_field = exchanged.support
                if sup == exchanged.support and dec == exchanged.dec:
                    sup_name = "nodes" if sup == MEDC.ON_NODES else "cells"
                    self.log(f"sending field {field_name!r} on {sup_name}...")
                    pfield = exchanged.para_field.getField()
                    array = field.getArray()
                    pfield.setArray(array)
                    pfield.setName(field_name)
                    self.log(repr(pfield), verbosity=2)
                    # self.log(array, verbosity=2)
                    dec.attachLocalField(pfield)
                    nb_field += 1
            if nb_field:
                self.log(f"sync...", verbosity=2)
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
        for sup, dec in self.dec.items():
            nb_field = 0
            for field_name in fields_names:
                exchanged = self.get_field(field_name)
                if sup == exchanged.support and dec == exchanged.dec:
                    sup_name = "nodes" if sup == MEDC.ON_NODES else "cells"
                    self.log(f"waiting for field {field_name!r} on {sup_name}...")
                    pfield = exchanged.para_field.getField()
                    pfield.setValues([0.0] * exchanged.size)
                    array = pfield.getArray()
                    array.rearrange(len(exchanged.components))
                    array.setInfoOnComponents(exchanged.components)
                    dec.attachLocalField(pfield)
                    self.log(repr(pfield), verbosity=2)
                    fields[field_name] = pfield
                    nb_field += 1
            if nb_field > 0:
                self.log(f"sync...", verbosity=2)
                dec.synchronize()
                self.log(f"recvData...", verbosity=2)
                dec.recvData()
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
        logger.debug("writing file {0!r}...".format(filename))
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
            *FieldOnNodesReal*: code_aster field defined on the whole mesh.
        """

        field = self._medcfield2aster(mc_field, field_type, time)
        return self.project_field(field)

    def export_field(self, field, field_name):
        """Convert a code_aster field defined on the whole mesh to
            a MEDCoupling field defined on the interface.

        Arguments:
            *FieldOnNodesReal*: code_aster field defined on the whole mesh.
            field_name (str): name of the field (like `DEPL`)

        Returns:
            *MEDCouplingField*: MEDCoupling field.
        """

        filename = "field_%i.med" % (MPI.ASTER_COMM_WORLD.rank)
        self._write_field2med(field, filename)

        fieldname = field.getName()[:8]
        mc_field = MEDC.ReadFieldNode(
            filename, field.getMesh().getName(), self.meshDimRelToMaxExt, fieldname, -1, -1
        )

        pfield = mc_field[self.interf_ids]
        pfield.setName(field_name)
        pfield.setNature(MEDC.IntensiveMaximum)
        pfield.setMesh(self.interf_mc)
        array = pfield.getArray()
        array.setInfoOnComponents(field.getComponents())
        pfield.checkConsistencyLight()
        os.remove(filename)

        return pfield

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
            displ (*FieldOnNodesReal*): code_aster displacement field.
            field_name (str): Field name.

        Returns:
            *MEDCouplingFieldDouble*: Displacement field.
        """
        cmps = [cmp for cmp in displ.getComponents() if cmp in ["DX", "DY", "DZ"]]
        if len(cmps) != len(displ.getComponents()):
            displ = displ.restrict(cmps=cmps)

        return self.export_field(displ, field_name)

    def export_temperature(self, temp, field_name):
        """Create a MEDCoupling field of temperature reduced on the interface mesh.

        Arguments:
            temp (*FieldOnNodesReal*): code_aster thermal field.
            field_name (str): Field name.

        Returns:
            *MEDCouplingFieldDouble*: Thermal field on cells.
        """

        cmps = [cmp for cmp in temp.getComponents() if cmp in ["TEMP"]]
        if len(cmps) != len(temp.getComponents()):
            temp = temp.restrict(cmps=cmps)

        return self.export_field(temp, field_name)

    def import_temperature(self, mc_temp, time=0.0):
        """Convert a MEDCoupling thermal field as a code_aster field.

        Arguments:
            mc_temp (*MEDCouplingFieldDouble*): MEDCoupling thermal field.
            time (float, optional): Time of assignment.

        Returns:
            *FieldOnNodesReal*: code_aster thermal field.
        """

        return self.import_field(mc_temp, "NOEU_TEMP_R", time)

    def export_pressure(self, pres, field_name):
        """Create a MEDCoupling field of pressure reduced on the interface mesh.

        Arguments:
            press (*FieldOnNodesReal*): code_aster pressure field.
            field_name (str): Field name.

        Returns:
            *MEDCouplingFieldDouble*: Pressure field on cells.
        """

        cmps = [cmp for cmp in pres.getComponents() if cmp in ["PRES"]]
        if len(cmps) != len(pres.getComponents()):
            pres = pres.restrict(cmps=cmps)

        return self.export_field(pres, field_name)

    def import_pressure(self, mc_pres, time=0.0):
        """Convert a MEDCoupling pressure field as a code_aster field.

        Arguments:
            mc_pres (*MEDCouplingFieldDouble*): MEDCoupling pressure field.
            time (float, optional): Time of assignment.

        Returns:
            *FieldOnNodesReal*: code_aster pressure field.
        """

        return self.import_field(mc_pres, "NOEU_PRES_R", time)

    def import_fluidforces(self, mc_fluidf, model, time=0.0):
        """Convert a MEDCoupling pressure field as a code_aster field.

        Arguments:
            mc_fluidf (*MEDCouplingField*): MEDCoupling pressure field.
            model (Model): Mechanical model.
            time (float, optional): Time of assignment.

        Returns:
            *LoadResult*: code_aster pressure as *LoadResult*.
        """

        if self.model_interf is None:
            modelization = "3D" if self.mesh_interf.getDimension() == 3 else "D_PLAN"

            self.model_interf = AFFE_MODELE(
                MAILLAGE=self.mesh_interf,
                AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION=modelization),
                DISTRIBUTION=_F(METHODE="CENTRALISE"),
            )

        mc_mesh = mc_fluidf.getMesh()
        tmpfile = "fort.78"
        MEDC.WriteUMesh(tmpfile, mc_mesh, True)
        MEDC.WriteFieldUsingAlreadyWrittenMesh(tmpfile, mc_fluidf)
        forc_elem = LIRE_CHAMP(
            MAILLAGE=self.mesh_interf,
            MODELE=self.model_interf,
            UNITE=78,
            NOM_MED=mc_fluidf.getName(),
            TYPE_CHAM="ELEM_FORC_R",
            NOM_CMP_IDEM="OUI",
        )
        os.remove(tmpfile)

        evol_char = LoadResult()
        evol_char.allocate(1)
        evol_char.setModel(self.model_interf, 0)
        evol_char.setTime(time, 0)
        evol_char.setField(forc_elem, "FORC_NODA", 0)

        if self.matr_proj is None:
            self.matr_proj = PROJ_CHAMP(
                METHODE="COLLOCATION", PROJECTION="NON", MODELE_1=self.model_interf, MODELE_2=model
            )

        return PROJ_CHAMP(
            RESULTAT=evol_char, MATR_PROJECTION=self.matr_proj, TOUT_ORDRE="OUI", MODELE_2=model
        )

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

        return PROJ_CHAMP(CHAM_GD=field, MATR_PROJECTION=self.matr_proj)
