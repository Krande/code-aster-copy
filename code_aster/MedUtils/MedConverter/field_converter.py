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

# person_in_charge: francesco.bettonte at edf.fr

import libaster
import numpy as np
from ...Utilities import medcoupling as medc, ParaMEDMEM as PMM


def getPhysicalQuantityFromFieldName(fname):
    """
    Get the aster physical quantity from the name of field.

    Physical Quantity has to be the final part of the field name.

    Arguments:
        name (str) : Field name.

    Returns :
        phys (str) : Physical Quantity name.

    """
    pqm = libaster.PhysicalQuantityManager.getAllPhysicalQuantityNames()
    candidates = [i for i in pqm if fname.endswith(i)]

    if not len(candidates) == 1:
        msg = f"Invalid field name (`{fname}`); expecting a physical quantity as end of name."
        raise RuntimeError(msg)
    else:
        name = candidates[0]
    return name


def getSymbolicNameFromMedField(medfield):
    """
    Get the aster symbolic name of field.

    Arguments:
        field (MEDFileField | MEDCouplingField) : The field in med format ( medcoupling ).

    Returns :
        phys (str) : Symbolic name of field.
        scal (str) : Scalar type of field.
    """
    trad = {medc.ON_NODES: "NOEU", medc.ON_CELLS: "ELEM", medc.ON_GAUSS_PT: "ELGA"}

    exceptions = {"DEPL_NOEU": "DEPL", "TEMP_NOEU": "TEMP", "PRES_NOEU": "PRES"}

    fname = medfield.getName()
    quantity, scal = getPhysicalQuantityFromFieldName(fname).split("_")

    phys = "_".join((quantity, trad[medfield.getTypeOfField()]))

    for src, dest in exceptions.items():
        phys = phys.replace(src, dest)

    return phys, scal


def toMCFieldAndProfileNode(asfield, medmesh, prefix=""):
    """Internal Function. Export the field to a new MEDCoupling field and profile.

    Arguments:
        asfield (*SimpleFieldOnNodes*): The aster field as simple field.
        medmesh (*MEDCouplingUMesh*): The medcoupling support mesh.
        prefix,  optional (str): Prefix for field names.

    Returns:
        field ( MEDCouplingFieldDouble ) : The field medcoupling format.

    """
    if not isinstance(medmesh, (medc.MEDFileUMesh, medc.MEDCouplingUMesh, PMM.MEDCouplingUMesh)):
        msg = "Argument must be a MEDCouplingUMesh, not '{}'"
        raise TypeError(msg.format(type(medmesh).__name__))

    # Aster values
    values, mask = asfield.toNumpy()

    # Restrict field based on mask
    restricted_nodes = np.where(np.any(mask, axis=1) == True)[0]
    restricted_values = values[restricted_nodes, :]

    # The field name containing the physical quantity
    field_name = "".join((prefix, asfield.getPhysicalQuantity()))
    profile_name = "".join((prefix, "NodesProfile"))

    # Medcoupling field
    field_values = medc.DataArrayDouble(restricted_values)
    field_values.setInfoOnComponents(asfield.getComponents())
    field_values.setName(field_name)

    medc_node_field = medc.MEDCouplingFieldDouble(medc.ON_NODES, medc.ONE_TIME)
    medc_node_field.setName(field_name)
    medc_node_field.setArray(field_values)
    medc_node_field.setNature(medc.IntensiveMaximum)

    # Med profile
    field_profile = medc.DataArrayInt(restricted_nodes)

    # Med support mesh for field ( restricted to profile nodes ) without cells
    field_mesh = medc.MEDCouplingUMesh()
    field_mesh.setName(medmesh.getName())
    field_mesh.setMeshDimension(medmesh.getMeshDimension())
    field_mesh.setCoords(medmesh.getCoords()[field_profile])
    field_mesh.allocateCells()
    medc_node_field.setMesh(field_mesh)

    field_profile.setName(profile_name)
    medc_node_field.checkConsistencyLight()

    return medc_node_field, field_profile


def toMCFieldAndProfileElem(asfield, medmesh, prefix=""):
    """Export the field to a new MEDCoupling field

    Arguments:
        medmesh (*MEDCouplingUMesh*): The medcoupling support mesh.

    Returns:
        field ( MEDCouplingFieldDouble ) : The field medcoupling format.
    """

    if not isinstance(medmesh, (medc.MEDCouplingUMesh, PMM.MEDCouplingUMesh)):
        msg = "toMedCouplingField() argument must be a MEDCouplingUMesh, not '{}'"
        raise TypeError(msg.format(type(medmesh).__name__))

    values, mask = asfield._cache["val"], asfield._cache["msk"]

    # Restrict field based on mask
    restricted_cells = np.where(np.any(mask, axis=1))[0]
    restricted_values = values[restricted_cells, :]

    # The field name containing the physical quantity
    field_name = "".join((prefix, asfield.getPhysicalQuantity()))
    profile_name = "".join((prefix, "ElemProfile"))

    # Medcoupling field
    field_values = medc.DataArrayDouble(restricted_values)
    field_values.setInfoOnComponents(asfield.getComponents())
    field_values.setName(field_name)
    medc_cell_field = medc.MEDCouplingFieldDouble(medc.ON_CELLS, medc.ONE_TIME)
    medc_cell_field.setName(field_name)
    medc_cell_field.setArray(field_values)
    medc_cell_field.setNature(medc.IntensiveConservation)

    field_profile = medc.DataArrayInt(restricted_cells)
    field_profile.setName("".join((prefix, "ElemProfile")))

    if len(restricted_cells) == medmesh.getNumberOfCells():
        medc_cell_field.setMesh(medmesh)
    else:
        raise NotImplementedError()

    medc_cell_field.checkConsistencyLight()

    return medc_cell_field, field_profile


def toMCFieldAndProfile(asfield, medmesh, prefix=""):
    loc = asfield.getLocalization()
    if loc == "NOEU":
        field, profile = toMCFieldAndProfileNode(asfield, medmesh, prefix)
    elif loc == "ELEM":
        field, profile = toMCFieldAndProfileElem(asfield, medmesh, prefix)
    else:
        raise NotImplementedError(loc)

    return field, profile


def toMedCouplingField(asfield, medmesh, prefix=""):
    """Export the field to a new MEDCoupling field

    Arguments:
        asfield (*SimpleFieldOnNodes*): The aster field as simple field.
        medmesh (*MEDCouplingUMesh*): The medcoupling support mesh.
        prefix,  optional (str): Prefix for field names.

    Returns:
        field ( MEDCouplingFieldDouble ) : The field medcoupling format.
    """

    medc_field, field_profile = toMCFieldAndProfile(asfield, medmesh, prefix)

    return medc_field


def toMedFileField1TS(asfield, medmesh, profile=False, prefix=""):
    """Export the field to a new MED field

    Arguments:
        asfield (*SimpleFieldOnNodes*): The aster field as simple field.
        medmesh (*MEDFileUMesh*): The medcoupling support mesh.
        profile, optional (bool): True to create a MED profile from field mask.
        prefix,  optional (str): Prefix for field names.

    Returns:
        field ( MEDFileField1TS ) : The field in med format ( medcoupling ).
    """
    medc_node_field, field_profile = toMCFieldAndProfile(asfield, medmesh, prefix)

    # Med field with profile
    medfield = medc.MEDFileField1TS()
    if profile:
        medfield.setFieldProfile(medc_node_field, medmesh, 1, field_profile)
    else:
        medfield.setFieldNoProfileSBT(medc_node_field)

    return medfield


def fromMedFileField1TSNodes(mc_field, astermesh):
    """Export the field information to create a new aster field.

    Arguments:
        field (*MEDCouplingFieldDouble*): The medcoupling field.
        mesh (Mesh) : field support mesh.

    Returns:
        phys (str) : Physical quantity name
        cmps (list[str]) : Name of components
        values (list[float]) : Field values

    """

    if not isinstance(mc_field, (medc.MEDCouplingFieldDouble, PMM.MEDCouplingFieldDouble)):
        msg = "Argument must be a MEDCouplingFieldDouble, not '{}'"
        raise TypeError(msg.format(type(mc_field).__name__))

    if mc_field.getTypeOfField() != medc.ON_NODES:
        raise RuntimeError("Field is not defined on nodes.")

    src, target = mc_field.getMesh().getNumberOfNodes(), astermesh.getNumberOfNodes()
    if src != target:
        raise RuntimeError(
            "Meshes have an incompatible number of nodes (%d != %d)." % (src, target)
        )

    arr = mc_field.getArray()
    cmps = arr.getInfoOnComponents()
    # FIXME : on devrait garder la physical quantity dans setDescription
    field_name = arr.getName() or mc_field.getName()
    phys = getPhysicalQuantityFromFieldName(field_name)
    values = arr.getValues()

    return phys, cmps, values


def fromMedFileField1TSCells(mc_field, astermesh):
    """Export the field information to create a new aster field.

    Arguments:
        field (*MEDCouplingFieldDouble*): The medcoupling field.
        mesh (Mesh) : field support mesh.

    Returns:
        phys (str) : Physical quantity name
        cmps (list[str]) : Name of components
        values (list[float]) : Field values

    """

    if not isinstance(mc_field, (medc.MEDCouplingFieldDouble, PMM.MEDCouplingFieldDouble)):
        msg = "Argument must be a MEDCouplingFieldDouble, not '{}'"
        raise TypeError(msg.format(type(mc_field).__name__))

    if mc_field.getTypeOfField() != medc.ON_CELLS:
        raise RuntimeError("Field is not defined on cells.")

    src, target = mc_field.getMesh().getNumberOfCells(), astermesh.getNumberOfCells()
    if src != target:
        raise RuntimeError(
            "Meshes have an incompatible number of cells (%d != %d)." % (src, target)
        )

    arr = mc_field.getArray()
    cmps = arr.getInfoOnComponents()
    # FIXME : on devrait garder la physical quantity dans setDescription
    field_name = arr.getName() or mc_field.getName()
    phys = getPhysicalQuantityFromFieldName(field_name)
    values = arr.getValues()

    return phys, cmps, values
