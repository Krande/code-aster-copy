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


def importMEDCDisplacement(mc_displ, modinterf, model, time=0.0):
    """Convert a MEDCoupling displacement field defined on the interface as
    a code_aster field.

    Arguments:
        mc_displ (*MEDCouplingField*): MEDCoupling displacement field.
        modinterf (*Model*): code_aster model of the interface.
        model (*Model*): code_aster model (entire).
        time (float, optional): Time of assignment.

    Returns:
        *FieldOnNodes*: code_aster displacement field.
    """
    mc_mesh = mc_displ.getMesh()
    tmpfile = "fort.77"
    MEDC.WriteUMesh(tmpfile, mc_mesh, True)
    MEDC.WriteFieldUsingAlreadyWrittenMesh(tmpfile, mc_displ)
    ca_displ = LIRE_CHAMP(
        UNITE=77,
        MAILLAGE=modinterf.getMesh(),
        # MODELE=modinterf,
        PROL_ZERO="OUI",
        NOM_MED="Displ",
        TYPE_CHAM="NOEU_DEPL_R",
        NOM_CMP_IDEM="OUI",
        INST=time,
        INFO=2,
    )
    os.remove(tmpfile)

    result = CREA_RESU(
        TYPE_RESU="EVOL_ELAS",
        OPERATION="AFFE",
        AFFE=_F(CHAM_GD=ca_displ, NOM_CHAM="DEPL", MODELE=modinterf, INST=time),
    )
    # projection on the entire model
    proj_result = PROJ_CHAMP(
        METHODE="COLLOCATION",
        RESULTAT=result,
        MODELE_1=modinterf,
        MODELE_2=model,
        # VIS_A_VIS=_F(
        #     GROUP_MA_1="interf",
        #     GROUP_MA_2="interf",
        # ),
        TOUT_ORDRE="OUI",
        PROL_ZERO="OUI",
    )
    displ = CREA_CHAMP(
        OPERATION="EXTR", RESULTAT=proj_result, NOM_CHAM="DEPL", TYPE_CHAM="NOEU_DEPL_R", INST=time
    )
    return displ


def exportMEDCPressure(temp, model, mesh_interf, ids):
    """Create a MEDCoupling field of pressure reduced on the interface mesh.

    Arguments:
        temp (*FieldOnNodes*): code_aster field.
        model (*Model*): full code_aster mechanical model (used for conversion).
        mesh_interf (*MEDCouplingUMesh*): mesh of the interface.
        ids (*DataArrayInt*): Numbering of interface cells.

    Returns:
        *MEDCouplingField*: Pressure field on cells.
    """
    nodal_press = CREA_CHAMP(
        OPERATION="ASSE",
        TYPE_CHAM="NOEU_PRES_R",
        MAILLAGE=model.getMesh(),
        ASSE=(_F(GROUP_MA="interf", CHAM_GD=temp, NOM_CMP="TEMP", NOM_CMP_RESU="PRES"),),
    )
    press = CREA_CHAMP(
        OPERATION="DISC",
        TYPE_CHAM="ELEM_PRES_R",
        MODELE=model,
        CHAM_GD=nodal_press,
        PROL_ZERO="OUI",
    )

    filename = "/tmp/press.med"
    _write_field2med(press, filename)

    mc_press = MEDC.ReadFieldCell(filename, press.getMesh().getName(), -1, press.getName(), -1, -1)
    pressure = mc_press[ids]
    pressure.setName("Pressure")
    pressure.setNature(MEDC.IntensiveConservation)
    pressure.setMesh(mesh_interf)
    pressure.getArray().setInfoOnComponents(["PRES [Pa]"])
    pressure.checkConsistencyLight()
    return pressure


def importMEDCPressure(mc_press, modinterf, model, time=0.0):
    """Convert a MEDCoupling pressure field as a code_aster field.

    Arguments:
        mc_press (*MEDCouplingField*): MEDCoupling pressure field.
        modinterf (*Model*): code_aster model of the interface.
        model (*Model*): full code_aster model.
        time (float, optional): Time of assignment.
        mesh (*Mesh*): code_aster mesh.

    Returns:
        *LoadResult*: code_aster pressure as *LoadResult*.
    """
    mc_mesh = mc_press.getMesh()
    tmpfile = "fort.78"
    MEDC.WriteUMesh(tmpfile, mc_mesh, True)
    MEDC.WriteFieldUsingAlreadyWrittenMesh(tmpfile, mc_press)
    press = LIRE_CHAMP(
        MAILLAGE=modinterf.getMesh(),
        MODELE=modinterf,
        UNITE=78,
        NOM_MED="Pressure",
        TYPE_CHAM="ELEM_PRES_R",
        NOM_CMP_IDEM="OUI",
        INST=time,
        INFO=2,
    )
    os.remove(tmpfile)

    press = CREA_RESU(
        TYPE_RESU="EVOL_CHAR",
        OPERATION="AFFE",
        AFFE=_F(NOM_CHAM="PRES", CHAM_GD=press, MODELE=modinterf, INST=time),
    )
    # projection on the entire model
    proj_press = PROJ_CHAMP(
        METHODE="COLLOCATION",
        RESULTAT=press,
        MODELE_1=modinterf,
        MODELE_2=model,
        # VIS_A_VIS=_F(
        #     GROUP_MA_1="interf",
        #     GROUP_MA_2="interf",
        # ),
        TOUT_ORDRE="OUI",
        PROL_ZERO="OUI",
    )
    return proj_press


def importMEDCFluidForces(mc_fluidf, field_name, modinterf, model, time=0.0):
    """Convert a MEDCoupling pressure field as a code_aster field.

    Arguments:
        mc_fluidf (*MEDCouplingField*): MEDCoupling pressure field.
        field_name (str): Field name.
        modinterf (*Model*): code_aster model of the interface.
        model (*Model*): full code_aster model.
        time (float, optional): Time of assignment.

    Returns:
        *LoadResult*: code_aster pressure as *LoadResult*.
    """
    mc_mesh = mc_fluidf.getMesh()
    tmpfile = "fort.78"
    MEDC.WriteUMesh(tmpfile, mc_mesh, True)
    MEDC.WriteFieldUsingAlreadyWrittenMesh(tmpfile, mc_fluidf)
    forces = LIRE_CHAMP(
        MAILLAGE=modinterf.getMesh(),
        MODELE=modinterf,
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
        AFFE=_F(NOM_CHAM="FORC_NODA", CHAM_GD=forces, MODELE=modinterf, INST=time),
    )
    # projection on the entire model
    proj_forces = PROJ_CHAMP(
        METHODE="COLLOCATION",
        RESULTAT=forces,
        MODELE_1=modinterf,
        MODELE_2=model,
        # VIS_A_VIS=_F(
        #     GROUP_MA_1="interf",
        #     GROUP_MA_2="interf",
        # ),
        TOUT_ORDRE="OUI",
        PROL_ZERO="OUI",
    )
    return proj_forces


def exportMEDCDisplacement(displ, field_name, mesh_interf, ids):
    """Create a MEDCoupling field of displacement reduced on the interface mesh.

    Arguments:
        displ (*FieldOnNodes*): code_aster displacement field.
        field_name (str): Field name.
        mesh_interf (*MEDCouplingUMesh*): mesh of the interface.
        ids (*DataArrayInt*): Numbering of interface cells.

    Returns:
        *MEDCouplingField*: Displacement field.
    """
    filename = "/tmp/displ.med"
    _write_field2med(displ, filename)

    fieldname = displ.getName()[:8]
    mc_displ = MEDC.ReadFieldNode(filename, displ.getMesh().getName(), -1, fieldname, -1, -1)

    displ = mc_displ[ids]
    displ.setName(field_name)
    displ.setNature(MEDC.IntensiveMaximum)
    displ.setMesh(mesh_interf)
    array = displ.getArray()
    array.setInfoOnComponents(["DX", "DY", "DZ"])
    displ.checkConsistencyLight()

    return displ


def _write_field2med(field, filename):
    if MPI.ASTER_COMM_WORLD.rank == 0:
        if os.path.exists(filename):
            os.remove(filename)
        print("writing file {0!r}...".format(filename), flush=True)
        field.printMedFile(filename)
    MPI.ASTER_COMM_WORLD.Barrier()
