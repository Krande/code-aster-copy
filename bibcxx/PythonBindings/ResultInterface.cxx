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

    MaterialFieldPtr ( Result::*c1 )() const =
        &Result::getMaterialField;
    MaterialFieldPtr ( Result::*c2 )( ASTERINTEGER ) const =
        &Result::getMaterialField;

    ModelPtr ( Result::*c3 )() const =
        &Result::getModel;
    ModelPtr ( Result::*c4 )( ASTERINTEGER ) const =
        &Result::getModel;

    ElementaryCharacteristicsPtr ( Result::*c5 )() const =
        &Result::getElementaryCharacteristics;
    ElementaryCharacteristicsPtr ( Result::*c6 )( ASTERINTEGER ) const =
        &Result::getElementaryCharacteristics;

    bool ( Result::*c7 )( const std::string ) const =
        &Result::printMedFile;
    bool ( Result::*c8 )( const std::string, std::string ) const =
        &Result::printMedFile;

    bool ( Result::*c9 )( const FieldOnNodesRealPtr,
                          const std::string&, const ASTERINTEGER ) = &Result::setField;
    bool ( Result::*c10 )( const FieldOnCellsRealPtr,
                          const std::string&, const ASTERINTEGER ) = &Result::setField;
    bool ( Result::*c21 )( const ConstantFieldOnCellsChar16Ptr,
                          const std::string&, const ASTERINTEGER ) = &Result::setField;

    void ( Result::*c11 )( const ModelPtr & ) = &Result::setModel;
    void ( Result::*c12 )( const ModelPtr &,  ASTERINTEGER ) = &Result::setModel;

    void ( Result::*c13 )( const ElementaryCharacteristicsPtr & ) =
        &Result::setElementaryCharacteristics;
    void ( Result::*c14 )( const ElementaryCharacteristicsPtr &,  ASTERINTEGER ) =
        &Result::setElementaryCharacteristics;

    void ( Result::*c15 )( const MaterialFieldPtr & ) =
        &Result::setMaterialField;
    void ( Result::*c16 )( const MaterialFieldPtr &,  ASTERINTEGER ) =
        &Result::setMaterialField;

    py::class_< Result, Result::ResultPtr,
            py::bases< DataStructure > >( "Result", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< Result, std::string >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< Result, std::string, std::string >))
        .def( "allocate", &Result::allocate, R"(
Allocate result

Arguments:
    nb_rank [int] :  number of rank to allocate

Returns:
    bool: True if allocation is ok
        )", ( py::arg("self" ), py::arg("nb_rank" )))
        .def( "setTimeValue", &Result::setTimeValue, R"(
Add time at the specified rank

Arguments:
    time [float] : time value to save
    rank [int] :  rank where to save time value
        )", ( py::arg("self" ), py::arg("time" ), py::arg("rank" )))
        .def( "getTimeValue", &Result::getTimeValue, R"(
Get time at the specified rank

Arguments:
    rank [int] :  rank where to save time value

Returns
    [float] : time value
        )", ( py::arg("self" ), py::arg("rank" )))
        .def( "addFieldOnNodesDescription", &Result::addFieldOnNodesDescription )
        .def( "setMaterialField", c15, R"(
Set material field on all ranks

Argument:
    mater [MaterialFieldPtr]: material field to set.
        )", ( py::arg("self" ), py::arg("mater")) )
        .def( "setMaterialField", c16, R"(
Set material field on the specified rank

Argument:
    mater [MaterialFieldPtr]: material field to set.
    rank [rank]: rank to set
        )", ( py::arg("self" ), py::arg("mater"), py::arg("rank")) )
        .def( "setListOfLoads", &Result::setListOfLoads, R"(
Set list of loads on the specified rank

Argument:
    load [ListOfLoadsPtr]: list of loads to set.
    rank [rank]: rank to set
        )", ( py::arg("self" ), py::arg("load"), py::arg("rank")) )
        .def( "setModel", c11, R"(
Set model on all ranks

Argument:
    model [ModelPtr]: model to set.
        )", ( py::arg("self" ), py::arg("model"))  )
        .def( "setModel", c12, R"(
Set model on the specified rank

Argument:
    model [ModelPtr]: model to set
    rank [rank]: rank to set
        )", ( py::arg("self" ), py::arg("model"), py::arg("rank"))  )
        .def( "setElementaryCharacteristics", c13,  R"(
Set elementary characterictics on all ranks

Argument:
    cara_elem [ElementaryCharacteristicsPtr]: elementary characterictics to set.
        )", ( py::arg("self" ), py::arg("mater")) )
        .def( "setElementaryCharacteristics", c14,  R"(
Set elementary characterictics on the specified rank

Argument:
    cara_elem [ElementaryCharacteristicsPtr]: elementary characterictics to set.
    rank [rank]: rank to set
        )", ( py::arg("self" ), py::arg("cara_elem"), py::arg("rank")) )
        .def( "printListOfFields", &Result::printListOfFields )
        .def( "getAllElementaryCharacteristics", &Result::getAllElementaryCharacteristics, R"(
Return the list of all elementary characteristics used in the result

Returns:
    list[ElementaryCharacteristicsPtr]: list of ElementaryCharacteristics.
        )", ( py::arg("self" )) )
        .def( "getElementaryCharacteristics", c5, R"(
Get elementary characterictics if only one is used else an execption is throw

Return:
    cara_elem [ElementaryCharacteristicsPtr]: a pointer to elementary characterictics.
        )", ( py::arg("self" )) )
        .def( "getElementaryCharacteristics", c6 , R"(
Get elementary characterictics at the specfied rank

Argument:
    rank [int]: rank

Return:
    cara_elem [ElementaryCharacteristicsPtr]: a pointer to elementary characterictics.
        )", ( py::arg("self" ), py::arg("rank")))
        .def( "getMaterialFields", &Result::getMaterialFields, R"(
Return the list of all material fields used in the result

Returns:
    list[MaterialFieldPtr]: list of material field.
        )", ( py::arg("self" )) )
        .def( "getMaterialField", c1 )
        .def( "getMaterialField", c2 )
        .def( "getListOfLoads", &Result::getListOfLoads, R"(
Get list of loads on the specified rank

Argument:
    rank [rank]: rank to get

Return:
    ListOfLoadsPtr: a pointer to list of loads.

        )", ( py::arg("self" ),  py::arg("rank")) )
        .def( "getMesh", &Result::getMesh, R"(
Return a pointer to mesh

Returns:
    mesh [BaseMeshPtr]: a pointer to the mesh.
        )", ( py::arg("self" )) )
        .def( "getModels", &Result::getModels, R"(
Return the list of all models used in the result

Returns:
    list[ModelPtr]: list of models.
        )", ( py::arg("self" )) )
        .def( "getModel", c3 )
        .def( "getModel", c4 )
        .def( "getNumberOfRanks", &Result::getNumberOfRanks )
        .def( "getAccessParameters", &Result::getAccessParameters, R"(
Return the access parameters of the result as Python dict.

Returns:
    dict{str : list[int,float,str]}: Dict of values for each access variable.
        )", ( py::arg("self" )))
        .def( "getFieldsOnNodesNames", &Result::getFieldsOnNodesNames, R"(
Return the names of the fields on nodes as Python list.

Returns:
    list[str]: List of names of the fields on nodes.
        )", ( py::arg("self" )))
        .def( "getFieldsOnCellsNames", &Result::getFieldsOnCellsNames, R"(
Return the names of the fields on cells as Python list.

Returns:
    list[str]: List of names of the fields on cells.
        )", ( py::arg("self" )))
        .def( "getConstantFieldsOnCellsNames", &Result::getConstantFieldsOnCellsNames, R"(
Return the names of the contant fields on cells as Python list.

Returns:
    list[str]: List of names of the contant fields on cells.
        )", ( py::arg("self" )))
        .def( "getRanks", &Result::getRanks, R"(
Return the list of ranks used to store fields

Returns:
    list[int]: List of ranks used to store fields.
        )", ( py::arg("self" )) )
        .def( "getFieldOnNodesReal", &Result::getFieldOnNodesReal )
        .def( "getFieldOnCellsReal", &Result::getFieldOnCellsReal )
        .def( "getConstantFieldOnCellsChar16", &Result::getConstantFieldOnCellsChar16 )
        .def( "printMedFile", c7 )
        .def( "printMedFile", c8 )
        .def( "setMesh", &Result::setMesh )
        .def( "build", &Result::build )
        .def( "printInfo", &Result::printInfo )
        .def( "setField", c9, R"(
Set a FieldOnNodes to result

Arguments:
    field [FieldOnNodesRealPtr] : field to set
    name [str]: symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    rank [int]: rank to set the field

Returns:
    bool: True if ok else False.
        )", ( py::arg("self"), py::arg("field"), py::arg("name"), py::arg("rank")) )
        .def( "setField", c10, R"(
Set a FieldOnCells to result

Arguments:
    field [FieldOnCellsRealPtr] : field to set
    name [str]: symbolic name of the field in the result (ex: 'VARI_ELGA', 'SIEF_ELGA'...)
    rank [int]: rank to set the field

Returns:
    bool: True if ok else False.
        )", ( py::arg("self"), py::arg("field"), py::arg("name"), py::arg("rank")) )
        .def( "setField", c21, R"(
Set a ConstantFieldOnCellsChar16 to result

Arguments:
    field [ConstantFieldOnCellsChar16] : field to set
    name [str]: symbolic name of the field in the result (ex: 'COMPOR', ...)
    rank [int]: rank to set the field

Returns:
    bool: True if ok else False.
        )", ( py::arg("self"), py::arg("field"), py::arg("name"), py::arg("rank")) )
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
