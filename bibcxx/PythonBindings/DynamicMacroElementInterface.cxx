/**
 * @file DynamicMacroElementInterface.cxx
 * @brief Interface python de DynamicMacroElement
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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

#include "PythonBindings/DynamicMacroElementInterface.h"

#include "aster_pybind.h"

void exportDynamicMacroElementToPython( py::module_ &mod ) {

    py::class_< DynamicMacroElement, DynamicMacroElement::DynamicMacroElementPtr, DataStructure >(
        mod, "DynamicMacroElement" )
        .def( py::init( &initFactoryPtr< DynamicMacroElement > ) )
        .def( py::init( &initFactoryPtr< DynamicMacroElement, std::string > ) )
        .def( "getDampingMatrix", &DynamicMacroElement::getDampingMatrix )
        .def( "getDOFNumbering", &DynamicMacroElement::getDOFNumbering )
        .def( "getImpedanceDampingMatrix", &DynamicMacroElement::getImpedanceDampingMatrix )
        .def( "getImpedanceMatrix", &DynamicMacroElement::getImpedanceMatrix )
        .def( "getImpedanceMassMatrix", &DynamicMacroElement::getImpedanceMassMatrix )
        .def( "getImpedanceStiffnessMatrix", &DynamicMacroElement::getImpedanceStiffnessMatrix )
        .def( "getMassMatrix", &DynamicMacroElement::getMassMatrix )
        .def( "getNumberOfNodes", &DynamicMacroElement::getNumberOfNodes )
        .def( "getStiffnessMatrixComplex", &DynamicMacroElement::getStiffnessMatrixComplex )
        .def( "getStiffnessMatrixReal", &DynamicMacroElement::getStiffnessMatrixReal )
        .def( "getMechanicalMode", &DynamicMacroElement::getMechanicalMode )
        .def( "getGeneralizedStiffnessMatrix", &DynamicMacroElement::getGeneralizedStiffnessMatrix )
        .def( "getGeneralizedMassMatrix", &DynamicMacroElement::getGeneralizedMassMatrix )
        .def( "getGeneralizedDampingMatrix", &DynamicMacroElement::getGeneralizedDampingMatrix )
        .def( "setDampingMatrix", &DynamicMacroElement::setDampingMatrix )
        .def( "setImpedanceDampingMatrix", &DynamicMacroElement::setImpedanceDampingMatrix )
        .def( "setImpedanceMatrix", &DynamicMacroElement::setImpedanceMatrix )
        .def( "setImpedanceMassMatrix", &DynamicMacroElement::setImpedanceMassMatrix )
        .def( "setImpedanceStiffnessMatrix", &DynamicMacroElement::setImpedanceStiffnessMatrix )
        .def( "setMassMatrix", &DynamicMacroElement::setMassMatrix )
        .def( "setMechanicalMode", &DynamicMacroElement::setMechanicalMode )
        .def( "setStiffnessMatrix",
              py::overload_cast< const AssemblyMatrixDisplacementComplexPtr & >(
                  &DynamicMacroElement::setStiffnessMatrix ) )
        .def( "setStiffnessMatrix", py::overload_cast< const AssemblyMatrixDisplacementRealPtr & >(
                                        &DynamicMacroElement::setStiffnessMatrix ) );
};
