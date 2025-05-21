/**
 * @file MedVectorInterface.cxx
 * @brief Interface python de MedVector
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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

#include "PythonBindings/MedVectorInterface.h"

#include "aster_pybind.h"

#ifdef ASTER_HAVE_MED
void exportMedVectorToPython( py::module_ &mod ) {

    py::class_< MedVector< double >, MedVectorPtr >( mod, "MedVector" )
        .def( py::init( &initFactoryPtr< MedVector< double > > ) )
        .def( py::init( &initFactoryPtr< MedVector< double >, int > ) )
        .def( define_pickling< MedVector< double > >() )

        .def( "getComponentName", &MedVector< double >::getComponentName, R"(
Get component name
            )" )
        .def( "getComponentNumber", &MedVector< double >::getComponentNumber, R"(
Get component name
            )" )
        .def( "getComponentVector", &MedVector< double >::getComponentVector, R"(
Get component on element vector
            )" )
        .def( "getCumulatedSizesVector", &MedVector< double >::getCumulatedSizesVector, R"(
Get cumulated sizes vector

Returns:
    list: Cumulated sizes for each element
            )" )
        .def( "getValues", &MedVector< double >::getValues, R"(
Get vector values (WARNING values are owned by MedVector: no copy)

Returns:
    numpy array: all field values
            )" )
        .def( "setComponentName", &MedVector< double >::setComponentName, R"(
Set component name
            )" )
        .def( "setComponentNumber", &MedVector< double >::setComponentNumber, R"(
Set component number
            )" )
        .def( "setComponentVector", &MedVector< double >::setComponentVector, R"(
Set component on element vector
            )" )
        .def( "setCumulatedSizesVector", &MedVector< double >::setCumulatedSizesVector, R"(
Set cumulated sizes vector
            )" )
        .def( "setValues", &MedVector< double >::setValues, R"(
Set vector values (WARNING values are owned by MedVector: no copy)
            )" )
        .def( "size", &MedVector< double >::size, R"(
Get vector size, ie: number of elements (cells or nodes)
)" );
}

#endif
