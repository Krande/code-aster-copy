/**
 * @file FieldOnNodesInterface.cxx
 * @brief Python interface for FieldOnNodes
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

/* person_in_charge: nicolas.sellenet at edf.fr */
#include "PythonBindings/FieldOnNodesInterface.h"

#include "aster_pybind.h"

#include "DataFields/MeshCoordinatesField.h"
#include "PythonBindings/DataStructureInterface.h"

void exportFieldOnNodesToPython( py::module_ &mod ) {
    /**
     * Object FieldOnNodesReal
     */

    py::class_< FieldOnNodesReal, FieldOnNodesRealPtr, DataField >( mod, "FieldOnNodesReal" )
        .def( py::init( &initFactoryPtr< FieldOnNodesReal > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesReal, std::string > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesReal, const FieldOnNodesReal & > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesReal, BaseDOFNumberingPtr > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesReal, MeshCoordinatesFieldPtr > ) )
        .def( "duplicate", &FieldOnNodesReal::duplicate )
        .def( "exportToSimpleFieldOnNodes", &FieldOnNodesReal::exportToSimpleFieldOnNodes )
        .def( "getPhysicalQuantity", &FieldOnNodesReal::getPhysicalQuantity )
        .def( "getMesh", &FieldOnNodesReal::getMesh )
        .def(
            "__getitem__",
            +[]( const FieldOnNodesReal &v, ASTERINTEGER i ) { return v.operator[]( i ); } )
        .def(
            "__setitem__",
            +[]( FieldOnNodesReal &v, ASTERINTEGER i, float f ) { return v.operator[]( i ) = f; } )
        .def( py::self += py::self )
        .def( py::self -= py::self )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( py::self * float() )
        .def( float() * py::self )
        .def( py::self *= float() )
        .def( -py::self )
        .def( "printMedFile", &FieldOnNodesReal::printMedFile, py::arg( "fileName" ),
              py::arg( "local" ) = true )
        .def( "setMesh", &FieldOnNodesReal::setMesh )
        .def( "setDescription", &FieldOnNodesReal::setDescription )
        .def( "build", &FieldOnNodesReal::build )
        .def( "getMesh", &FieldOnNodesReal::getMesh )
        .def( "getDescription", &FieldOnNodesReal::getDescription )
        .def( "updateValuePointers", &FieldOnNodesReal::updateValuePointers )
        .def( "getComponents", &FieldOnNodesReal::getComponents, R"(
            Get list of components

            Returns:
                list[str]: list of components
            )" )
        .def( "getNumberOfComponents", &FieldOnNodesReal::getNumberOfComponents, R"(
            Get number of components

            Returns:
                int: number of components
            )" )
        .def( "norm", &FieldOnNodesReal::norm, R"(
            Return the euclidean norm of the field

            Arguments:
                normType (str): "NORM_1", "NORM_2", "NORM_INFINITY"
                list_cmp (list[str]) : list of components used to compute norm (default: all)

            Returns:
                float: euclidean norm
            )",
              py::arg( "normType" ), py::arg( "list_cmp" ) = VectorString() )
        .def( "dot", &FieldOnNodesReal::dot, R"(
            Return the dot product of two fields

            Arguments:
                other (FieldOnNodes): other field

            Returns:
                float: dot product
            )",
              py::arg( "other" ) )
        .def( "size", &FieldOnNodesReal::size, R"(
            Return the size of the field

            Returns:
                int: number of element in the field
            )" )
        .def( "scale", &FieldOnNodesReal::scale, R"(
            Scale in-place the field by a diagonal matrix stored as an array

            Arguments:
                vect (float): diagonal matrix stored as an array
            )",
              py::arg( "vect" ) )
        .def( "setValues", py::overload_cast< const ASTERDOUBLE & >( &FieldOnNodesReal::setValues ),
              R"(
            Set values of the field

            Arguments:
                value (float): value to set
            )",
              py::arg( "value" ) )
        .def( "setValues", py::overload_cast< const VectorReal & >( &FieldOnNodesReal::setValues ),
              R"(
            Set values of the field

            Arguments:
                values (list[float]): list of values to set
            )",
              py::arg( "values" ) )
        .def( "setValues",
              py::overload_cast< const std::map< std::string, ASTERDOUBLE > & >(
                  &FieldOnNodesReal::setValues ),
              R"(
            Set values of the field where components and values are given as a dict.
            If the component is not present in the field then it is discarded
            Example: { "X1" : 0.0, "X3" : 0.0 }

            Arguments:
                value (dict[str, float]): dict of values to set (key: str, value: float)
            )",
              py::arg( "value" ) )
        .def( "getValues", &FieldOnNodesReal::getValues, py::return_value_policy::reference, R"(
            Return a list of values as (x1, y1, z1, x2, y2, z2...)

            Returns:
                list[float]: List of values.
            )" )
        .def( "getNodesAndComponentsNumberFromDOF",
              &FieldOnNodesReal::getNodesAndComponentsNumberFromDOF, R"(
            Return a list of values such that for each DOF, it gives the node id and component id
            as [dof1=[node_1, comp_1], dof2=[node_1, comp_2], ....]

            Returns:
                list[[int, int]]: List of values (node, component) for each DOF.
            )" );
    /**
     * Object FieldOnNodesComplex
     */
    py::class_< FieldOnNodesComplex, FieldOnNodesComplexPtr, DataField >( mod,
                                                                          "FieldOnNodesComplex" )
        .def( py::init( &initFactoryPtr< FieldOnNodesComplex > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesComplex, std::string > ) )
        .def( py::init< const FieldOnNodesComplex & >() )
        .def( py::init( &initFactoryPtr< FieldOnNodesComplex, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesComplex, BaseDOFNumberingPtr > ) )
        .def( "exportToSimpleFieldOnNodes", &FieldOnNodesComplex::exportToSimpleFieldOnNodes )
        .def( "getPhysicalQuantity", &FieldOnNodesComplex::getPhysicalQuantity )
        .def( "getMesh", &FieldOnNodesComplex::getMesh )
        .def(
            "__getitem__",
            +[]( const FieldOnNodesComplex &v, ASTERINTEGER i ) { return v.operator[]( i ); } )
        .def(
            "__setitem__", +[]( FieldOnNodesComplex &v, ASTERINTEGER i,
                                ASTERCOMPLEX f ) { return v.operator[]( i ) = f; } )
        .def( "printMedFile", &FieldOnNodesComplex::printMedFile )
        .def( "setMesh", &FieldOnNodesComplex::setMesh )
        .def( "setDescription", &FieldOnNodesComplex::setDescription )
        .def( "build", &FieldOnNodesComplex::build )
        .def( "getMesh", &FieldOnNodesComplex::getMesh )
        .def( "getDescription", &FieldOnNodesComplex::getDescription )
        .def( "getValues", &FieldOnNodesComplex::getValues, py::return_value_policy::reference, R"(
            Return a list of complex values as [x11, x21, ..., xm1, x12, x22, ..., xm2...]
            (m is the total number of componenets)

            Returns:
                list[complex]: List of values.
            )" )
        .def( "getComponents", &FieldOnNodesComplex::getComponents, R"(
            Get list of components

            Returns:
                list[str]: list of components
            )" )
        .def( "getNumberOfComponents", &FieldOnNodesComplex::getNumberOfComponents, R"(
            Get number of components

            Returns:
                int: number of components
            )" )
        .def( "scale", &FieldOnNodesComplex::scale, R"(
            Scale in-place the field by a diagonal matrix stored as an array

            Arguments:
                vect (float): diagonal matrix stored as an array
            )",
              py::arg( "vect" ) )
        .def( "dot", &FieldOnNodesComplex::dot, R"(
            Return the dot product of two complex fields

            Arguments:
                other (FieldOnNodes): other field

            Returns:
                complex: dot product
            )",
              py::arg( "other" ) )
        .def( "norm", &FieldOnNodesComplex::norm, R"(
            Return the euclidean norm of the field

            Arguments:
                normType (str): "NORM_1", "NORM_2", "NORM_INFINITY"
                list_cmp (list[str]) : list of components used to compute norm (default: all)

            Returns:
                float: euclidean norm
            )",
              py::arg( "normType" ), py::arg( "list_cmp" ) = VectorString() )
        .def( "setValues",
              py::overload_cast< const ASTERCOMPLEX & >( &FieldOnNodesComplex::setValues ), R"(
            Set values of the field

            Arguments:
                value (complex): value to set
            )",
              py::arg( "value" ) )
        .def( "setValues",
              py::overload_cast< const VectorComplex & >( &FieldOnNodesComplex::setValues ),
              R"(
            Set values of the field

            Arguments:
                values (list[complex]): list of values to set
            )",
              py::arg( "values" ) )
        .def( "updateValuePointers", &FieldOnNodesComplex::updateValuePointers );

    /**
     * Object FieldOnNodesLong
     */
    py::class_< FieldOnNodesLong, FieldOnNodesLongPtr, DataField >( mod, "FieldOnNodesLong" )
        .def( py::init( &initFactoryPtr< FieldOnNodesLong > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesLong, std::string > ) )
        .def( py::init< const FieldOnNodesLong & >() )
        .def( py::init( &initFactoryPtr< FieldOnNodesLong, BaseDOFNumberingPtr > ) )
        .def( "setDescription", &FieldOnNodesLong::setDescription )
        .def( "getDescription", &FieldOnNodesLong::getDescription )
        .def( "getMesh", &FieldOnNodesLong::getMesh )
        .def( "setMesh", &FieldOnNodesLong::setMesh );

    /**
     * Object FieldOnNodesChar8
     */
    py::class_< FieldOnNodesChar8, FieldOnNodesChar8Ptr, DataField >( mod, "FieldOnNodesChar8" )
        .def( py::init( &initFactoryPtr< FieldOnNodesChar8 > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesChar8, std::string > ) )
        .def( py::init< const FieldOnNodesChar8 & >() )
        .def( py::init( &initFactoryPtr< FieldOnNodesChar8, BaseDOFNumberingPtr > ) )
        .def( "setDescription", &FieldOnNodesChar8::setDescription )
        .def( "getDescription", &FieldOnNodesChar8::getDescription )
        .def( "getMesh", &FieldOnNodesChar8::getMesh )
        .def( "setMesh", &FieldOnNodesChar8::setMesh );
};
