/**
 * @file ResultInterface.cxx
 * @brief Interface python de Result
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
 *
 *   This file is part of Code_Aster.
 *
 *   Code_Aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Code_Aster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "PythonBindings/ResultInterface.h"

#include "aster_pybind.h"

void exportResultToPython( py::module_ &mod ) {

    py::class_< Result, Result::ResultPtr, DataStructure >( mod, "Result" )
        .def( py::init( &initFactoryPtr< Result, std::string > ) )
        .def( py::init( &initFactoryPtr< Result, std::string, std::string > ) )
        .def( "allocate", &Result::allocate, R"(
Allocate result

Arguments:
    nb_rank (int):  number of rank to allocate
        )",
              py::arg( "nb_rank" ) )
        .def( "setTimeValue", &Result::setTimeValue, R"(
Add time at the specified rank

Arguments:
    time (float): time value to save
    rank (int):  rank where to save time value
        )",
              py::arg( "time" ), py::arg( "rank" ) )
        .def( "setParameterValue", &Result::setParameterValue, R"(
Add theta at the specified rank

Arguments:
    name (float): parameter name to store
    value (float): parameter value to store
    rank (int):  rank where to save time value
        )",
              py::arg( "para_name" ), py::arg( "value" ), py::arg( "rank" ) )
        .def( "getTimeValue", &Result::getTimeValue, R"(
Get time at the specified rank

Arguments:
    rank (int):  rank where to save time value

Returns
    float: time value
        )",
              py::arg( "rank" ) )
        .def( "addFieldOnNodesDescription", &Result::addFieldOnNodesDescription )
        .def( "setMaterialField",
              py::overload_cast< const MaterialFieldPtr & >( &Result::setMaterialField ), R"(
Set material field on all ranks

Arguments:
    mater (MaterialField): material field to set.
        )",
              py::arg( "mater" ) )
        .def( "setMaterialField",
              py::overload_cast< const MaterialFieldPtr &, ASTERINTEGER >(
                  &Result::setMaterialField ),
              R"(
Set material field on the specified rank

Arguments:
    mater (MaterialField): material field to set.
    rank (int): rank to set
        )",
              py::arg( "mater" ), py::arg( "rank" ) )
        .def( "setListOfLoads", &Result::setListOfLoads, R"(
Set list of loads on the specified rank

Arguments:
    load (ListOfLoads): list of loads to set.
    rank (int): rank to set
        )",
              py::arg( "load" ), py::arg( "rank" ) )
        .def( "setModel", py::overload_cast< const ModelPtr & >( &Result::setModel ), R"(
Set model on all ranks

Arguments:
    model (Model): model to set.
        )",
              py::arg( "model" ) )
        .def( "setModel", py::overload_cast< const ModelPtr &, ASTERINTEGER >( &Result::setModel ),
              R"(
Set model on the specified rank

Arguments:
    model (Model): model to set
    rank (int): rank to set
        )",
              py::arg( "model" ), py::arg( "rank" ) )
        .def( "setElementaryCharacteristics",
              py::overload_cast< const ElementaryCharacteristicsPtr & >(
                  &Result::setElementaryCharacteristics ),
              R"(
Set elementary characterictics on all ranks

Arguments:
    cara_elem (ElementaryCharacteristics): elementary characterictics to set.
        )",
              py::arg( "cara_elem" ) )
        .def( "setElementaryCharacteristics",
              py::overload_cast< const ElementaryCharacteristicsPtr &, ASTERINTEGER >(
                  &Result::setElementaryCharacteristics ),
              R"(
Set elementary characterictics on the specified rank

Arguments:
    cara_elem (ElementaryCharacteristics): elementary characterictics to set.
    rank (int): rank to set
        )",
              py::arg( "cara_elem" ), py::arg( "rank" ) )
        .def( "printListOfFields", &Result::printListOfFields, R"(
Print the names of all fields (real, complex, ...) stored in the result.
        )" )
        .def( "getAllElementaryCharacteristics", &Result::getAllElementaryCharacteristics, R"(
Return the list of all elementary characteristics used in the result

Returns:
    list[ElementaryCharacteristics]: list of ElementaryCharacteristics.
        )" )
        .def(
            "getElementaryCharacteristics",
            py::overload_cast< ASTERINTEGER >( &Result::getElementaryCharacteristics, py::const_ ),
            R"(
Get elementary characterictics at the specfied rank

Arguments:
    rank (int): rank

Returns:
    ElementaryCharacteristics: a pointer to elementary characterictics.
        )",
            py::arg( "rank" ) )
        .def( "getElementaryCharacteristics",
              py::overload_cast<>( &Result::getElementaryCharacteristics, py::const_ ) )
        .def(
            "hasElementaryCharacteristics",
            py::overload_cast< ASTERINTEGER >( &Result::hasElementaryCharacteristics, py::const_ ),
            R"(
Test if a elementary characterictics is used at the specfied rank

Arguments:
    rank (int): rank

Returns:
    bool: *True* if at least one elementary characterictics used else *False*.
        )",
            py::arg( "rank" ) )
        .def( "hasElementaryCharacteristics",
              py::overload_cast<>( &Result::hasElementaryCharacteristics, py::const_ ) )
        .def( "hasListOfLoads",
              py::overload_cast< const ASTERINTEGER & >( &Result::hasListOfLoads, py::const_ ), R"(
Test if a list of loads is used at the specfied rank

Arguments:
    rank (int): rank

Returns:
    bool: *True* if at least one list of loads is used else *False*.
        )",
              py::arg( "rank" ) )
        .def( "hasListOfLoads", py::overload_cast<>( &Result::hasListOfLoads, py::const_ ) )
        .def( "hasModel", &Result::hasModel, R"(
Test if a model is used at the specfied rank

Arguments:
    rank (int): rank

Returns:
    bool: *True* if at a model used else *False*.
        )",
              py::arg( "rank" ) )
        .def( "hasMaterialField", &Result::hasMaterialField, R"(
Test if a material field is used at the specfied rank

Arguments:
    rank (int): rank

Returns:
    bool: *True* if at a material field used else *False*.
        )",
              py::arg( "rank" ) )
        .def( "getMaterialFields", &Result::getMaterialFields, R"(
Return the list of all material fields used in the result

Returns:
    list[MaterialField]: list of material field.
        )" )
        .def( "getMaterialField",
              py::overload_cast< ASTERINTEGER >( &Result::getMaterialField, py::const_ ),
              R"(
Return the material field for the given rank.

Arguments:
    rank (int): rank

Returns:
    MaterialField: Material field.
              )",
              py::arg( "rank" ) )
        .def( "getMaterialField", py::overload_cast<>( &Result::getMaterialField, py::const_ ) )
        .def( "getListOfLoads", &Result::getListOfLoads, R"(
Get list of loads on the specified rank

Arguments:
    rank (int): rank to get

Returns:
    ListOfLoads: a pointer to list of loads.

        )",
              py::arg( "rank" ) )
        .def( "getMesh", &Result::getMesh, R"(
Return a pointer to mesh

Returns:
    mesh (Mesh): a pointer to the mesh.
        )" )
        .def( "getModels", &Result::getModels, R"(
Return the list of all models used in the result

Returns:
    list[Model]: list of models.
        )" )
        .def( "getModel", py::overload_cast< ASTERINTEGER >( &Result::getModel, py::const_ ), R"(
Return the model for the given rank.

Arguments:
    rank (int): rank

Returns:
    Model: Model object.
              )",
              py::arg( "rank" ) )
        .def( "getModel", py::overload_cast<>( &Result::getModel, py::const_ ) )
        .def( "getNumberOfRanks", &Result::getNumberOfRanks, R"(
Get the number of rank stored in the result

Returns:
    int: number of rank stored.
        )" )
        .def( "getAccessParameters", &Result::getAccessParameters, R"(
Return the access parameters of the result as Python dict.

Returns:
    dict{str : list[int,float,str]}: Dict of values for each access variable.
        )" )
        .def( "getFieldsOnNodesRealNames", &Result::getFieldsOnNodesRealNames, R"(
Return the names of the real fields on nodes as Python list.

Returns:
    list(str): List of names of the real fields on nodes.
        )" )
        .def( "getFieldsOnCellsRealNames", &Result::getFieldsOnCellsRealNames, R"(
Return the names of the real fields on cells as Python list.

Returns:
    list(str): List of names of the real fields on cells.
        )" )
        .def( "getFieldsOnNodesComplexNames", &Result::getFieldsOnNodesComplexNames, R"(
Return the names of the complex fields on nodes as Python list.

Returns:
    list(str): List of names of the complex fields on nodes.
        )" )
        .def( "getFieldsOnCellsComplexNames", &Result::getFieldsOnCellsComplexNames, R"(
Return the names of the complex fields on cells as Python list.

Returns:
    list(str): List of names of the complex fields on cells.
        )" )
        .def( "getFieldsOnCellsLongNames", &Result::getFieldsOnCellsLongNames, R"(
Return the names of the integer fields on cells as Python list.

Returns:
    list(str): List of names of the integer fields on cells.
        )" )
        .def( "getConstantFieldsOnCellsChar16Names", &Result::getConstantFieldsOnCellsChar16Names,
              R"(
Return the names of the contant char16 fields on cells as Python list.

Returns:
    list(str): List of names of the contant fields on cells.
        )" )
        .def( "getConstantFieldsOnCellsRealNames", &Result::getConstantFieldsOnCellsRealNames,
              R"(
Return the names of the contant real fields on cells as Python list.

Returns:
    list(str): List of names of the contant fields on cells.
        )" )
        .def( "getRanks", &Result::getRanks, R"(
Return the list of ranks used to store fields

Returns:
    list[int]: List of ranks used to store fields.
        )" )
        .def( "getFieldOnNodesReal", &Result::getFieldOnNodesReal, R"(
Get a FieldOnNodesReal from result.

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    rank (int): rank to set the field

Returns:
    FieldOnNodesReal: field to get
        )",
              py::arg( "name" ), py::arg( "rank" ) )
        .def( "getFieldOnCellsReal", &Result::getFieldOnCellsReal, R"(
Get a FieldOnCellsReal from result.

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    rank (int): rank to set the field

Returns:
    FieldOnCellsReal: field to get
        )",
              py::arg( "name" ), py::arg( "rank" ) )
        .def( "getFieldOnNodesComplex", &Result::getFieldOnNodesComplex, R"(
Get a FieldOnNodesComplex from result.

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    rank (int): rank to set the field

Returns:
    FieldOnNodesComplex: field to get
        )",
              py::arg( "name" ), py::arg( "rank" ) )
        .def( "getFieldOnCellsComplex", &Result::getFieldOnCellsComplex, R"(
Get a FieldOnCellsComplex from result.

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    rank (int): rank to set the field

Returns:
    FieldOnCellsComplex: field to get
        )",
              py::arg( "name" ), py::arg( "rank" ) )
        .def( "getFieldOnCellsLong", &Result::getFieldOnCellsLong, R"(
Get a FieldOnCellsLong from result.

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    rank (int): rank to set the field

Returns:
    FieldOnCellsLong: field to get
        )",
              py::arg( "name" ), py::arg( "rank" ) )
        .def( "getConstantFieldOnCellsChar16", &Result::getConstantFieldOnCellsChar16, R"(
Get a ConstantFieldOnCellsChar16 from result.

Arguments:
    name (str): symbolic name of the field in the result (ex: 'COMPORTEMENT', ...)
    rank (int): rank to set the field

Returns:
    ConstantFieldOnCellsChar16: field to get
        )",
              py::arg( "name" ), py::arg( "rank" ) )
        .def( "getConstantFieldOnCellsReal", &Result::getConstantFieldOnCellsReal, R"(
Get a ConstantFieldOnCellsReal from result.

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    rank (int): rank to set the field

Returns:
    ConstantFieldOnCellsReal: field to get
        )",
              py::arg( "name" ), py::arg( "rank" ) )
        .def( "printMedFile", &Result::printMedFile,
              R"(
Print the result in a MED file.

Args:
    filename (str): Path to the output file.
    medname (str): Name of the result in the MED file. (default: "")
    local (bool): Print only the local domain if *True*. (default: True)
              )",
              py::arg( "filename" ), py::arg( "medname" ) = "", py::arg( "local" ) = true )
        .def( "setMesh", &Result::setMesh, R"(
Set the mesh used by the result.

Arguments:
    mesh (BaseMesh): mesh to set

        )",
              py::arg( "mesh" ) )
        .def( "build", &Result::build, R"(
Build the result from the name of the result. It stores fields which are setted in c++ or
created in fortran

Arguments:
    feds (list[FiniteElementDescriptor]) : list of additional finite element descriptor used to
        build FieldOnCells
    fnds (list[FieldOnNodesDescriptionPtr]) : list of additional field description used to
        build FieldOnNodes

Returns:
    bool: *True* if ok.
        )",
              py::arg( "feds" ) = std::vector< FiniteElementDescriptorPtr >(),
              py::arg( "fnds" ) = std::vector< FieldOnNodesDescriptionPtr >() )
        .def( "printInfo", &Result::printInfo )
        .def( "getFieldsNames", &Result::getFieldsNames, R"(
Return the list of names of stored fields

Returns:
    list[str]: List of names of stored fields.
        )" )
        .def( "resize", &Result::resize, R"(
Resize the object.

Arguments:
    nbRanks (int): new expected size. Should be greater than the current size,
        otherwise the size is unchanged.
        )",
              py::arg( "nbRanks" ) )
        .def(
            "setField",
            py::overload_cast< const FieldOnNodesRealPtr, const std::string &, const ASTERINTEGER >(
                &Result::setField ),
            R"(
Set a real FieldOnNodes to result.

Arguments:
    field (FieldOnNodesReal): field to set
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    rank (int): rank to set the field
        )",
            py::arg( "field" ), py::arg( "name" ), py::arg( "rank" ) )
        .def( "setField",
              py::overload_cast< const FieldOnNodesComplexPtr, const std::string &,
                                 const ASTERINTEGER >( &Result::setField ),
              R"(
Set a complex FieldOnNodes to result.

Arguments:
    field (FieldOnNodesComplex): field to set
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    rank (int): rank to set the field
        )",
              py::arg( "field" ), py::arg( "name" ), py::arg( "rank" ) )
        .def(
            "setField",
            py::overload_cast< const FieldOnCellsRealPtr, const std::string &, const ASTERINTEGER >(
                &Result::setField ),
            R"(
Set a real FieldOnCells to result

Arguments:
    field (FieldOnCellsReal): field to set
    name (str): symbolic name of the field in the result (ex: 'VARI_ELGA', 'SIEF_ELGA'...)
    rank (int): rank to set the field
        )",
            py::arg( "field" ), py::arg( "name" ), py::arg( "rank" ) )
        .def( "setField",
              py::overload_cast< const FieldOnCellsComplexPtr, const std::string &,
                                 const ASTERINTEGER >( &Result::setField ),
              R"(
Set a complex FieldOnCells to result

Arguments:
    field (FieldOnCellsComplex): field to set
    name (str): symbolic name of the field in the result (ex: 'VARI_ELGA', 'SIEF_ELGA'...)
    rank (int): rank to set the field
        )",
              py::arg( "field" ), py::arg( "name" ), py::arg( "rank" ) )
        .def(
            "setField",
            py::overload_cast< const FieldOnCellsLongPtr, const std::string &, const ASTERINTEGER >(
                &Result::setField ),
            R"(
Set a long FieldOnCells to result

Arguments:
    field (FieldOnCellsLong): field to set
    name (str): symbolic name of the field in the result (ex: 'VARI_ELGA', 'SIEF_ELGA'...)
    rank (int): rank to set the field
        )",
            py::arg( "field" ), py::arg( "name" ), py::arg( "rank" ) )
        .def( "setField",
              py::overload_cast< const ConstantFieldOnCellsChar16Ptr, const std::string &,
                                 const ASTERINTEGER >( &Result::setField ),
              R"(
Set a ConstantFieldOnCellsChar16 to result

Arguments:
    field (ConstantFieldOnCellsChar16): field to set
    name (str): symbolic name of the field in the result (ex: 'COMPOR', ...)
    rank (int): rank to set the field
        )",
              py::arg( "field" ), py::arg( "name" ), py::arg( "rank" ) )
        .def( "setField",
              py::overload_cast< const ConstantFieldOnCellsRealPtr, const std::string &,
                                 const ASTERINTEGER >( &Result::setField ),
              R"(
Set a ConstantFieldOnCellsReal to result

Arguments:
    field (ConstantFieldOnCellsReal): field to set
    name (str): symbolic name of the field in the result (ex: 'COMPOR', ...)
    rank (int): rank to set the field
        )",
              py::arg( "field" ), py::arg( "name" ), py::arg( "rank" ) )
        .def( "getTable", &ListOfTables::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              py::arg( "identifier" ) )
        .def( "getFieldsNames", &Result::getFieldsNames, R"(
Return the list of names of stored fields

Returns:
    list[str]: List of names of stored fields.
        )" )
        .def( "getFiniteElementDescriptors", &Result::getFiniteElementDescriptors, R"(
Get list of finite element descriptor to build internal FieldOnCells

Returns:
    list[FiniteElementDescriptor]: list of finite element descriptor
        )" )
        .def( "getFieldOnNodesDescriptions", &Result::getFieldOnNodesDescriptions, R"(
Get list of field's description to build internal FieldOnNodes

Returns:
    list[FieldOnNodesDescription]: list of field's description
        )" );
};
