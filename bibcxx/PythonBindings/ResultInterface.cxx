/**
 * @file ResultInterface.cxx
 * @brief Interface python de Result
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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
#include "PythonBindings/factory.h"
#include <boost/python.hpp>

namespace py = boost::python;

void exportResultToPython() {

    MaterialFieldPtr ( Result::*c1 )() const = &Result::getMaterialField;
    MaterialFieldPtr ( Result::*c2 )( ASTERINTEGER ) const = &Result::getMaterialField;

    ModelPtr ( Result::*c3 )() const = &Result::getModel;
    ModelPtr ( Result::*c4 )( ASTERINTEGER ) const = &Result::getModel;

    ElementaryCharacteristicsPtr ( Result::*c5 )() const = &Result::getElementaryCharacteristics;
    ElementaryCharacteristicsPtr ( Result::*c6 )( ASTERINTEGER ) const =
        &Result::getElementaryCharacteristics;

    bool ( Result::*c7 )( const std::string, bool ) const = &Result::printMedFile;
    bool ( Result::*c8 )( const std::string, std::string, bool ) const = &Result::printMedFile;
    bool ( Result::*c71 )( const std::string ) const = &Result::printMedFile;
    bool ( Result::*c81 )( const std::string, std::string ) const = &Result::printMedFile;

    bool ( Result::*c9 )( const FieldOnNodesRealPtr, const std::string &, const ASTERINTEGER ) =
        &Result::setField;
    bool ( Result::*c19 )( const FieldOnNodesComplexPtr, const std::string &, const ASTERINTEGER ) =
        &Result::setField;

    bool ( Result::*c10 )( const FieldOnCellsRealPtr, const std::string &, const ASTERINTEGER ) =
        &Result::setField;
    bool ( Result::*c20 )( const FieldOnCellsComplexPtr, const std::string &, const ASTERINTEGER ) =
        &Result::setField;
    bool ( Result::*c21 )( const FieldOnCellsLongPtr, const std::string &, const ASTERINTEGER ) =
        &Result::setField;

    bool ( Result::*c22 )( const ConstantFieldOnCellsChar16Ptr, const std::string &,
                           const ASTERINTEGER ) = &Result::setField;
    bool ( Result::*c23 )( const ConstantFieldOnCellsRealPtr, const std::string &,
                           const ASTERINTEGER ) = &Result::setField;

    void ( Result::*c11 )( const ModelPtr & ) = &Result::setModel;
    void ( Result::*c12 )( const ModelPtr &, ASTERINTEGER ) = &Result::setModel;

    void ( Result::*c13 )( const ElementaryCharacteristicsPtr & ) =
        &Result::setElementaryCharacteristics;
    void ( Result::*c14 )( const ElementaryCharacteristicsPtr &, ASTERINTEGER ) =
        &Result::setElementaryCharacteristics;

    void ( Result::*c15 )( const MaterialFieldPtr & ) = &Result::setMaterialField;
    void ( Result::*c16 )( const MaterialFieldPtr &, ASTERINTEGER ) = &Result::setMaterialField;

    bool ( Result::*c17 )() const = &Result::hasElementaryCharacteristics;
    bool ( Result::*c18 )( ASTERINTEGER ) const = &Result::hasElementaryCharacteristics;

    py::class_< Result, Result::ResultPtr, py::bases< DataStructure > >( "Result", py::no_init )
        .def( "__init__", py::make_constructor( &initFactoryPtr< Result, std::string > ) )
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< Result, std::string, std::string > ) )
        .def( "allocate", &Result::allocate, R"(
Allocate result

Arguments:
    nb_rank (int):  number of rank to allocate

Returns:
    bool: True if allocation is ok
        )",
              ( py::arg( "self" ), py::arg( "nb_rank" ) ) )
        .def( "setTimeValue", &Result::setTimeValue, R"(
Add time at the specified rank

Arguments:
    time (float): time value to save
    rank (int):  rank where to save time value
        )",
              ( py::arg( "self" ), py::arg( "time" ), py::arg( "rank" ) ) )
        .def( "getTimeValue", &Result::getTimeValue, R"(
Get time at the specified rank

Arguments:
    rank (int):  rank where to save time value

Returns
    float: time value
        )",
              ( py::arg( "self" ), py::arg( "rank" ) ) )
        .def( "addFieldOnNodesDescription", &Result::addFieldOnNodesDescription )
        .def( "setMaterialField", c15, R"(
Set material field on all ranks

Arguments:
    mater (MaterialField): material field to set.
        )",
              ( py::arg( "self" ), py::arg( "mater" ) ) )
        .def( "setMaterialField", c16, R"(
Set material field on the specified rank

Arguments:
    mater (MaterialField): material field to set.
    rank (int): rank to set
        )",
              ( py::arg( "self" ), py::arg( "mater" ), py::arg( "rank" ) ) )
        .def( "setListOfLoads", &Result::setListOfLoads, R"(
Set list of loads on the specified rank

Arguments:
    load (ListOfLoads): list of loads to set.
    rank (int): rank to set
        )",
              ( py::arg( "self" ), py::arg( "load" ), py::arg( "rank" ) ) )
        .def( "setModel", c11, R"(
Set model on all ranks

Arguments:
    model (Model): model to set.
        )",
              ( py::arg( "self" ), py::arg( "model" ) ) )
        .def( "setModel", c12, R"(
Set model on the specified rank

Arguments:
    model (Model): model to set
    rank (int): rank to set
        )",
              ( py::arg( "self" ), py::arg( "model" ), py::arg( "rank" ) ) )
        .def( "setElementaryCharacteristics", c13, R"(
Set elementary characterictics on all ranks

Arguments:
    cara_elem (ElementaryCharacteristics): elementary characterictics to set.
        )",
              ( py::arg( "self" ), py::arg( "mater" ) ) )
        .def( "setElementaryCharacteristics", c14, R"(
Set elementary characterictics on the specified rank

Arguments:
    cara_elem (ElementaryCharacteristics): elementary characterictics to set.
    rank (int): rank to set
        )",
              ( py::arg( "self" ), py::arg( "cara_elem" ), py::arg( "rank" ) ) )
        .def( "printListOfFields", &Result::printListOfFields, R"(
Print the names of all fields (real, complex, ...) stored in the result.
        )",
              ( py::arg( "self" ) ) )
        .def( "getAllElementaryCharacteristics", &Result::getAllElementaryCharacteristics, R"(
Return the list of all elementary characteristics used in the result

Returns:
    list[ElementaryCharacteristics]: list of ElementaryCharacteristics.
        )",
              ( py::arg( "self" ) ) )
        .def( "getElementaryCharacteristics", c5, R"(
Get elementary characterictics if only one is used else an execption is throw

Returns:
    ElementaryCharacteristics: a pointer to elementary characterictics.
        )",
              ( py::arg( "self" ) ) )
        .def( "getElementaryCharacteristics", c6, R"(
Get elementary characterictics at the specfied rank

Arguments:
    rank (int): rank

Returns:
    ElementaryCharacteristics: a pointer to elementary characterictics.
        )",
              ( py::arg( "self" ), py::arg( "rank" ) ) )
        .def( "hasElementaryCharacteristics", c17, R"(
Test if at least one elementary characterictics used

Returns:
    bool: *True* if at least one elementary characterictics used else *False*.
        )",
              ( py::arg( "self" ) ) )
        .def( "hasElementaryCharacteristics", c18, R"(
Test if a elementary characterictics is used at the specfied rank

Arguments:
    rank (int): rank

Returns:
    bool: *True* if at a elementary characterictics used else *False*.
        )",
              ( py::arg( "self" ), py::arg( "rank" ) ) )
        .def( "hasModel", &Result::hasModel, R"(
Test if a model is used at the specfied rank

Arguments:
    rank (int): rank

Returns:
    bool: *True* if at a model used else *False*.
        )",
              ( py::arg( "self" ), py::arg( "rank" ) ) )
        .def( "hasMaterialField", &Result::hasMaterialField, R"(
Test if a material field is used at the specfied rank

Arguments:
    rank (int): rank

Returns:
    bool: *True* if at a material field used else *False*.
        )",
              ( py::arg( "self" ), py::arg( "rank" ) ) )
        .def( "getMaterialFields", &Result::getMaterialFields, R"(
Return the list of all material fields used in the result

Returns:
    list[MaterialField]: list of material field.
        )",
              ( py::arg( "self" ) ) )
        .def( "getMaterialField", c1 )
        .def( "getMaterialField", c2 )
        .def( "getListOfLoads", &Result::getListOfLoads, R"(
Get list of loads on the specified rank

Arguments:
    rank (int): rank to get

Returns:
    ListOfLoads: a pointer to list of loads.

        )",
              ( py::arg( "self" ), py::arg( "rank" ) ) )
        .def( "getMesh", &Result::getMesh, R"(
Return a pointer to mesh

Returns:
    mesh (Mesh): a pointer to the mesh.
        )",
              ( py::arg( "self" ) ) )
        .def( "getModels", &Result::getModels, R"(
Return the list of all models used in the result

Returns:
    list[Model]: list of models.
        )",
              ( py::arg( "self" ) ) )
        .def( "getModel", c3 )
        .def( "getModel", c4 )
        .def( "getNumberOfRanks", &Result::getNumberOfRanks, R"(
Get the number of rank stored in the result

Returns:
    int: number of rank stored.
        )",
              ( py::arg( "self" ) ) )
        .def( "getAccessParameters", &Result::getAccessParameters, R"(
Return the access parameters of the result as Python dict.

Returns:
    dict{str : list[int,float,str]}: Dict of values for each access variable.
        )",
              ( py::arg( "self" ) ) )
        .def( "getFieldsOnNodesRealNames", &Result::getFieldsOnNodesRealNames, R"(
Return the names of the real fields on nodes as Python list.

Returns:
    list(str): List of names of the real fields on nodes.
        )",
              ( py::arg( "self" ) ) )
        .def( "getFieldsOnCellsRealNames", &Result::getFieldsOnCellsRealNames, R"(
Return the names of the real fields on cells as Python list.

Returns:
    list(str): List of names of the real fields on cells.
        )",
              ( py::arg( "self" ) ) )
        .def( "getFieldsOnNodesComplexNames", &Result::getFieldsOnNodesComplexNames, R"(
Return the names of the complex fields on nodes as Python list.

Returns:
    list(str): List of names of the complex fields on nodes.
        )",
              ( py::arg( "self" ) ) )
        .def( "getFieldsOnCellsComplexNames", &Result::getFieldsOnCellsComplexNames, R"(
Return the names of the complex fields on cells as Python list.

Returns:
    list(str): List of names of the complex fields on cells.
        )",
              ( py::arg( "self" ) ) )
        .def( "getFieldsOnCellsLongNames", &Result::getFieldsOnCellsLongNames, R"(
Return the names of the integer fields on cells as Python list.

Returns:
    list(str): List of names of the integer fields on cells.
        )",
              ( py::arg( "self" ) ) )
        .def( "getConstantFieldsOnCellsChar16Names", &Result::getConstantFieldsOnCellsChar16Names,
              R"(
Return the names of the contant char16 fields on cells as Python list.

Returns:
    list(str): List of names of the contant fields on cells.
        )",
              ( py::arg( "self" ) ) )
        .def( "getConstantFieldsOnCellsRealNames", &Result::getConstantFieldsOnCellsRealNames,
              R"(
Return the names of the contant real fields on cells as Python list.

Returns:
    list(str): List of names of the contant fields on cells.
        )",
              ( py::arg( "self" ) ) )
        .def( "getRanks", &Result::getRanks, R"(
Return the list of ranks used to store fields

Returns:
    list[int]: List of ranks used to store fields.
        )",
              ( py::arg( "self" ) ) )
        .def( "getFieldOnNodesReal", &Result::getFieldOnNodesReal, R"(
Get a FieldOnNodesReal from result.

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    rank (int): rank to set the field

Returns:
    FieldOnNodesReal: field to get
        )",
              ( py::arg( "self" ), py::arg( "name" ), py::arg( "rank" ) ) )
        .def( "getFieldOnCellsReal", &Result::getFieldOnCellsReal, R"(
Get a FieldOnCellsReal from result.

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    rank (int): rank to set the field

Returns:
    FieldOnCellsReal: field to get
        )",
              ( py::arg( "self" ), py::arg( "name" ), py::arg( "rank" ) ) )
        .def( "getFieldOnNodesComplex", &Result::getFieldOnNodesComplex, R"(
Get a FieldOnNodesComplex from result.

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    rank (int): rank to set the field

Returns:
    FieldOnNodesComplex: field to get
        )",
              ( py::arg( "self" ), py::arg( "name" ), py::arg( "rank" ) ) )
        .def( "getFieldOnCellsComplex", &Result::getFieldOnCellsComplex, R"(
Get a FieldOnCellsComplex from result.

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    rank (int): rank to set the field

Returns:
    FieldOnCellsComplex: field to get
        )",
              ( py::arg( "self" ), py::arg( "name" ), py::arg( "rank" ) ) )
        .def( "getFieldOnCellsLong", &Result::getFieldOnCellsLong, R"(
Get a FieldOnCellsLong from result.

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    rank (int): rank to set the field

Returns:
    FieldOnCellsLong: field to get
        )",
              ( py::arg( "self" ), py::arg( "name" ), py::arg( "rank" ) ) )
        .def( "getConstantFieldOnCellsChar16", &Result::getConstantFieldOnCellsChar16, R"(
Get a ConstantFieldOnCellsChar16 from result.

Arguments:
    name (str): symbolic name of the field in the result (ex: 'COMPORTEMENT', ...)
    rank (int): rank to set the field

Returns:
    ConstantFieldOnCellsChar16: field to get
        )",
              ( py::arg( "self" ), py::arg( "name" ), py::arg( "rank" ) ) )
        .def( "getConstantFieldOnCellsReal", &Result::getConstantFieldOnCellsReal, R"(
Get a ConstantFieldOnCellsReal from result.

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    rank (int): rank to set the field

Returns:
    ConstantFieldOnCellsReal: field to get
        )",
              ( py::arg( "self" ), py::arg( "name" ), py::arg( "rank" ) ) )
        .def( "printMedFile", c7 )
        .def( "printMedFile", c8 )
        .def( "printMedFile", c71 )
        .def( "printMedFile", c81 )
        .def( "setMesh", &Result::setMesh, R"(
Set the mesh used by the result.

Arguments:
    mesh (BaseMesh): mesh to set

        )",
              ( py::arg( "self" ), py::arg( "mesh" ) ) )
        .def( "build", &Result::build, R"(
Build the result from the name of the result. It stores fields which are setted in c++ or
created in fortran

Returns:
    bool: List of names of stored fields.
        )",
              ( py::arg( "self" ) ) )
        .def( "printInfo", &Result::printInfo )
        .def( "getFieldsNames", &Result::getFieldsNames, R"(
Return the list of names of stored fields

Returns:
    list[str]: List of names of stored fields.
        )",
              ( py::arg( "self" ) ) )
        .def( "resize", &Result::resize, R"(
Resize the object.

Arguments:
    nbRanks (int): new expected size. Should be greater than the current size,
        otherwise the size is unchanged.

Returns:
    bool: True if ok else False.
        )",
              ( py::arg( "self" ), py::arg( "nbRanks" ) ) )
        .def( "setField", c9, R"(
Set a real FieldOnNodes to result.

Arguments:
    field (FieldOnNodesReal): field to set
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    rank (int): rank to set the field

Returns:
    bool: True if ok else False.
        )",
              ( py::arg( "self" ), py::arg( "field" ), py::arg( "name" ), py::arg( "rank" ) ) )
        .def( "setField", c19, R"(
Set a complex FieldOnNodes to result.

Arguments:
    field (FieldOnNodesComplex): field to set
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    rank (int): rank to set the field

Returns:
    bool: True if ok else False.
        )",
              ( py::arg( "self" ), py::arg( "field" ), py::arg( "name" ), py::arg( "rank" ) ) )
        .def( "setField", c10, R"(
Set a real FieldOnCells to result

Arguments:
    field (FieldOnCellsReal): field to set
    name (str): symbolic name of the field in the result (ex: 'VARI_ELGA', 'SIEF_ELGA'...)
    rank (int): rank to set the field

Returns:
    bool: True if ok else False.
        )",
              ( py::arg( "self" ), py::arg( "field" ), py::arg( "name" ), py::arg( "rank" ) ) )
        .def( "setField", c20, R"(
Set a complex FieldOnCells to result

Arguments:
    field (FieldOnCellsComplex): field to set
    name (str): symbolic name of the field in the result (ex: 'VARI_ELGA', 'SIEF_ELGA'...)
    rank (int): rank to set the field

Returns:
    bool: True if ok else False.
        )",
              ( py::arg( "self" ), py::arg( "field" ), py::arg( "name" ), py::arg( "rank" ) ) )
        .def( "setField", c21, R"(
Set a long FieldOnCells to result

Arguments:
    field (FieldOnCellsLong): field to set
    name (str): symbolic name of the field in the result (ex: 'VARI_ELGA', 'SIEF_ELGA'...)
    rank (int): rank to set the field

Returns:
    bool: True if ok else False.
        )",
              ( py::arg( "self" ), py::arg( "field" ), py::arg( "name" ), py::arg( "rank" ) ) )
        .def( "setField", c22, R"(
Set a ConstantFieldOnCellsChar16 to result

Arguments:
    field (ConstantFieldOnCellsChar16): field to set
    name (str): symbolic name of the field in the result (ex: 'COMPOR', ...)
    rank (int): rank to set the field

Returns:
    bool: True if ok else False.
        )",
              ( py::arg( "self" ), py::arg( "field" ), py::arg( "name" ), py::arg( "rank" ) ) )
        .def( "setField", c23, R"(
Set a ConstantFieldOnCellsReal to result

Arguments:
    field (ConstantFieldOnCellsReal): field to set
    name (str): symbolic name of the field in the result (ex: 'COMPOR', ...)
    rank (int): rank to set the field

Returns:
    bool: True if ok else False.
        )",
              ( py::arg( "self" ), py::arg( "field" ), py::arg( "name" ), py::arg( "rank" ) ) )
        .def( "getTable", &ListOfTables::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              ( py::arg( "self" ), py::arg( "identifier" ) ) )

        ;
};
