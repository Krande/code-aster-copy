/**
 * @file StudyDescriptionInterface.cxx
 * @brief Interface python de StudyDescription
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

#include <boost/python.hpp>

namespace py = boost::python;
#include "PythonBindings/factory.h"
#include "PythonBindings/StudyDescriptionInterface.h"
#include "PythonBindings/LoadUtilities.h"

void exportStudyDescriptionToPython() {

    py::class_< StudyDescriptionClass, StudyDescriptionPtr > c1( "StudyDescription",
                                                                    py::no_init );
    // fake initFactoryPtr: not a DataStructure
    c1.def( "__init__",
            py::make_constructor(
                &initFactoryPtr< StudyDescriptionClass, ModelPtr, MaterialFieldPtr >));
    c1.def( "__init__",
            py::make_constructor(
                &initFactoryPtr< StudyDescriptionClass, ModelPtr, MaterialFieldPtr,
                                ElementaryCharacteristicsPtr >));
    c1.def( "getModel", &StudyDescriptionClass::getModel,
        py::return_value_policy<py::copy_const_reference>(), R"(
Return the model

Returns:
    ModelPtr: a pointer to the model
        )",
              ( py::arg( "self" ) )   );
    c1.def( "getMaterialField", &StudyDescriptionClass::getMaterialField,
        py::return_value_policy<py::copy_const_reference>(), R"(
Return the material field

Returns:
    MaterialFieldPtr: a pointer to the material field
        )",
              ( py::arg( "self" ) )   );
    c1.def( "getCodedMaterial", &StudyDescriptionClass::getCodedMaterial,
        py::return_value_policy<py::copy_const_reference>(), R"(
Return the coded material

Returns:
    CodedMaterialPtr: a pointer to the coded material
        )",
              ( py::arg( "self" ) )   );
    c1.def( "getElementaryCharacteristics",
        &StudyDescriptionClass::getElementaryCharacteristics,
        py::return_value_policy<py::copy_const_reference>(), R"(
Return the elementary charateristics

Returns:
    ElementaryCharacteristicsPtr: a pointer to the elementary charateristics
        )",
              ( py::arg( "self" ) )   );
    addDirichletBCToInterface( c1 );
    addMechanicalLoadToInterface( c1 );
#ifdef ASTER_HAVE_MPI
    addParallelMechanicalLoadToInterface( c1 );
#endif /* ASTER_HAVE_MPI */
};
