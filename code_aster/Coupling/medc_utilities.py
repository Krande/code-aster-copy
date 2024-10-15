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
Definition of utilities to convert fields.
"""

import os

import medcoupling as MEDC
import numpy as np

from ..Commands import LIRE_CHAMP, PROJ_CHAMP, CREA_RESU, CREA_CHAMP
from ..Utilities.mpi_utils import MPI


class MEDProj:
    """Class handling the MEDCoupling import/export.

    Arguments:
        mesh_interf (*MEDCouplingUMesh*): mesh of the interface.
        ids (*DataArrayInt*): Numbering of interface cells.
        meshDimRelToMaxExt (int) : a relative dimension of the mesh interface
            compared to whole mesh.
        modin_terf (*Model*): code_aster model of the interface.
        model (*Model*): code_aster model (entire).
    """

    def __init__(self, mesh_interf, ids, meshDimRelToMaxExt, model_interf, model):
        self.mesh_interf = mesh_interf
        self.ids = ids
        self.meshDimRelToMaxExt = meshDimRelToMaxExt
        self.model_interf = model_interf
        self.model = model

    @staticmethod
    def _write_field2med(field, filename):
        """Write field on disk using MED format.

        Arguments:
            field (FieldOn*): aster field.
            filename (str): name of MED file.
        """
        if MPI.ASTER_COMM_WORLD.rank == 0:
            if os.path.exists(filename):
                os.remove(filename)
            print("writing file {0!r}...".format(filename), flush=True)
            field.printMedFile(filename)
        MPI.ASTER_COMM_WORLD.Barrier()

    def importMEDCDisplacement(self, mc_displ, time=0.0):
        """Convert a MEDCoupling displacement field defined on the interface as
        a code_aster field.

        Arguments:
            mc_displ (*MEDCouplingField*): MEDCoupling displacement field.
            time (float, optional): Time of assignment.

        Returns:
            *FieldOnNodes*: code_aster displacement field.
        """
        mc_mesh = mc_displ.getMesh()
        tmpfile = "fort.77"
        MEDC.WriteUMesh(tmpfile, mc_mesh, True)
        MEDC.WriteFieldUsingAlreadyWrittenMesh(tmpfile, mc_displ)
        depl = LIRE_CHAMP(
            UNITE=77,
            MAILLAGE=self.model_interf.getMesh(),
            PROL_ZERO="OUI",
            NOM_MED=mc_displ.getName(),
            TYPE_CHAM="NOEU_DEPL_R",
            NOM_CMP_IDEM="OUI",
            INST=time,
            INFO=1,
        )
        os.remove(tmpfile)

        # projection on the entire model
        proj_depl = PROJ_CHAMP(
            METHODE="COLLOCATION",
            CHAM_GD=depl,
            MODELE_1=self.model_interf,
            MODELE_2=self.model,
            PROL_ZERO="OUI",
        )

        return proj_depl

    def exportMEDCDisplacement(self, displ, field_name):
        """Create a MEDCoupling field of displacement reduced on the interface mesh.

        Arguments:
            displ (*FieldOnNodes*): code_aster displacement field.
            field_name (str): Field name.

        Returns:
            *MEDCouplingField*: Displacement field.
        """
        filename = "/tmp/displ.med"
        self._write_field2med(displ, filename)
        cmps = [cmp for cmp in displ.getComponents() if cmp in ["DX", "DY", "DZ"]]

        fieldname = displ.getName()[:8]
        mc_displ = MEDC.ReadFieldNode(
            filename, displ.getMesh().getName(), self.meshDimRelToMaxExt, fieldname, -1, -1
        )

        displ = mc_displ[self.ids]
        displ.setName(field_name)
        displ.setNature(MEDC.IntensiveMaximum)
        displ.setMesh(self.mesh_interf)
        array = displ.getArray()
        array.setInfoOnComponents(cmps)
        displ.checkConsistencyLight()
        os.remove(filename)

        return displ

    def exportMEDCTemperature(self, temp, field_name):
        """Create a MEDCoupling field of temperature reduced on the interface mesh.

        Arguments:
            temp (*FieldOnNodes*): code_aster thermal field.
            field_name (str): Field name.

        Returns:
            *MEDCouplingField*: Pressure field on cells.
        """

        filename = "/tmp/temp.med"
        self._write_field2med(temp, filename)

        fieldname = temp.getName()[:8]
        mc_temp = MEDC.ReadFieldNode(
            filename, temp.getMesh().getName(), self.meshDimRelToMaxExt, fieldname, -1, -1
        )

        temp = mc_temp[self.ids]
        temp.setName(field_name)
        temp.setNature(MEDC.IntensiveMaximum)
        temp.setMesh(self.mesh_interf)
        temp.getArray().setInfoOnComponents(["TEMP"])
        temp.checkConsistencyLight()
        os.remove(filename)

        return temp

    def importMEDCTemperature(self, mc_temp, time=0.0):
        """Convert a MEDCoupling pemperature field as a code_aster field.

        Arguments:
            mc_temp (*MEDCouplingField*): MEDCoupling temp field.
            time (float, optional): Time of assignment.

        Returns:
            *FieldOnNodes*: code_aster thermal field.
        """
        mc_mesh = mc_temp.getMesh()
        tmpfile = "fort.78"
        MEDC.WriteUMesh(tmpfile, mc_mesh, True)
        MEDC.WriteFieldUsingAlreadyWrittenMesh(tmpfile, mc_temp)
        temp = LIRE_CHAMP(
            MAILLAGE=self.model_interf.getMesh(),
            PROL_ZERO="OUI",
            UNITE=78,
            NOM_MED=mc_temp.getName(),
            TYPE_CHAM="NOEU_TEMP_R",
            NOM_CMP_IDEM="OUI",
            INST=time,
            INFO=1,
        )
        os.remove(tmpfile)

        # projection on the whole model
        proj_temp = PROJ_CHAMP(
            METHODE="COLLOCATION",
            CHAM_GD=temp,
            MODELE_1=self.model_interf,
            MODELE_2=self.model,
            PROL_ZERO="OUI",
        )
        return proj_temp

    def importMEDCFluidForces(self, mc_fluidf, field_name, time=0.0):
        """Convert a MEDCoupling pressure field as a code_aster field.

        Arguments:
            mc_fluidf (*MEDCouplingField*): MEDCoupling pressure field.
            field_name (str): Field name.
            time (float, optional): Time of assignment.

        Returns:
            *LoadResult*: code_aster pressure as *LoadResult*.
        """
        mc_mesh = mc_fluidf.getMesh()
        tmpfile = "fort.78"
        MEDC.WriteUMesh(tmpfile, mc_mesh, True)
        MEDC.WriteFieldUsingAlreadyWrittenMesh(tmpfile, mc_fluidf)
        forces = LIRE_CHAMP(
            MAILLAGE=self.model_interf.getMesh(),
            MODELE=self.model_interf,
            UNITE=78,
            NOM_MED=field_name,
            TYPE_CHAM="ELEM_FORC_R",
            NOM_CMP_IDEM="OUI",
            INST=time,
            INFO=2,
        )
        os.remove(tmpfile)

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


def createField(mesh, name, components, fieldtype, val=0.0):
    """Create a field on nodes, filled by zero.

    Arguments:
       mesh (*MEDCouplingUMesh*): Mesh support.
        name (str): Field name.
        components (list[str]): List of components.
        val (float): value for initialization (default: 0.)

    Returns:
        *MEDCouplingField*: field.
    """
    field = MEDC.MEDCouplingFieldDouble(fieldtype, MEDC.ONE_TIME)
    field.setMesh(mesh)
    field.setName(name)
    field.setTime(0.0, 0, 0)
    field.setNature(
        MEDC.IntensiveMaximum if fieldtype == MEDC.ON_NODES else MEDC.IntensiveConservation
    )
    array = MEDC.DataArrayDouble()
    if fieldtype == MEDC.ON_NODES:
        size = mesh.getNumberOfNodes()
    elif fieldtype == MEDC.ON_CELLS:
        size = mesh.getNumberOfCells()
    else:
        raise RuntimeError("Unknown type.")
    array.alloc(size, len(components))
    array.setInfoOnComponents(components)
    array.fillWithValue(val)
    field.setArray(array)
    field.checkConsistencyLight()
    return field


def mergeComponents(name, *by_comp):
    """Merge fields containing one component into a field with all components.

    Arguments:
        name (str): Field name.
        by_comp (list[*MEDCouplingField*]): List of one component fields.

    Returns:
        *MEDCouplingField*: Field containing all components.
    """
    first = by_comp[0]
    mesh = first.getMesh()
    components, values = [], []
    for i in by_comp:
        arr = i.getArray()
        components.extend(arr.getInfoOnComponents())
        values.append(arr.toNumPyArray())
    field = MEDC.MEDCouplingFieldDouble(first.getTypeOfField(), MEDC.ONE_TIME)
    field.setName(name)
    field.setMesh(mesh)
    field.setTime(*first.getTime())
    field.setNature(first.getNature())
    array = MEDC.DataArrayDouble(np.concatenate(values))
    array.rearrange(len(components))
    array.setInfoOnComponents(components)
    field.setArray(array)
    field.checkConsistencyLight()
    return field
