/**
 * @file ExternalStateVariablesInterface.cxx
 * @brief Interface python de BaseExternalStateVariables
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

#include "PythonBindings/ListOfExternalStateVariablesInterface.h"

#include "aster_pybind.h"

#include "Materials/ListOfExternalStateVariables.h"

void exportListOfExternalStateVariablesToPython( py::module_ &mod ) {

    py::class_< ListOfExternalStateVariables, ListOfExternalStateVariablesPtr > c3(
        mod, "ListOfExternalStateVariables" );
    c3.def( py::init( &initFactoryPtr< ListOfExternalStateVariables, const MeshPtr & > ) );
    c3.def( py::init( &initFactoryPtr< ListOfExternalStateVariables, const SkeletonPtr & > ) );
#ifdef ASTER_HAVE_MPI
    c3.def( py::init( &initFactoryPtr< ListOfExternalStateVariables, const ParallelMeshPtr & > ) );
#endif /* ASTER_HAVE_MPI */
    c3.def( "addExternalStateVariableOnMesh",
            &ListOfExternalStateVariables::addExternalStateVariableOnMesh<
                TemperatureExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnGroupOfCells",
            &ListOfExternalStateVariables::addExternalStateVariableOnGroupOfCells<
                TemperatureExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnMesh",
            &ListOfExternalStateVariables::addExternalStateVariableOnMesh<
                GeometryExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnGroupOfCells",
            &ListOfExternalStateVariables::addExternalStateVariableOnGroupOfCells<
                GeometryExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnMesh",
            &ListOfExternalStateVariables::addExternalStateVariableOnMesh<
                CorrosionExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnGroupOfCells",
            &ListOfExternalStateVariables::addExternalStateVariableOnGroupOfCells<
                CorrosionExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnMesh",
            &ListOfExternalStateVariables::addExternalStateVariableOnMesh<
                IrreversibleDeformationExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnGroupOfCells",
            &ListOfExternalStateVariables::addExternalStateVariableOnGroupOfCells<
                IrreversibleDeformationExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnMesh",
            &ListOfExternalStateVariables::addExternalStateVariableOnMesh<
                ConcreteHydratationExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnGroupOfCells",
            &ListOfExternalStateVariables::addExternalStateVariableOnGroupOfCells<
                ConcreteHydratationExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnMesh",
            &ListOfExternalStateVariables::addExternalStateVariableOnMesh<
                IrradiationExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnGroupOfCells",
            &ListOfExternalStateVariables::addExternalStateVariableOnGroupOfCells<
                IrradiationExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnMesh",
            &ListOfExternalStateVariables::addExternalStateVariableOnMesh<
                SteelPhasesExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnGroupOfCells",
            &ListOfExternalStateVariables::addExternalStateVariableOnGroupOfCells<
                SteelPhasesExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnMesh",
            &ListOfExternalStateVariables::addExternalStateVariableOnMesh<
                ZircaloyPhasesExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnGroupOfCells",
            &ListOfExternalStateVariables::addExternalStateVariableOnGroupOfCells<
                ZircaloyPhasesExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnMesh",
            &ListOfExternalStateVariables::addExternalStateVariableOnMesh<
                Neutral1ExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnGroupOfCells",
            &ListOfExternalStateVariables::addExternalStateVariableOnGroupOfCells<
                Neutral1ExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnMesh",
            &ListOfExternalStateVariables::addExternalStateVariableOnMesh<
                Neutral2ExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnGroupOfCells",
            &ListOfExternalStateVariables::addExternalStateVariableOnGroupOfCells<
                Neutral2ExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnMesh",
            &ListOfExternalStateVariables::addExternalStateVariableOnMesh<
                Neutral3ExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnGroupOfCells",
            &ListOfExternalStateVariables::addExternalStateVariableOnGroupOfCells<
                Neutral3ExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnMesh",
            &ListOfExternalStateVariables::addExternalStateVariableOnMesh<
                ConcreteDryingExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnGroupOfCells",
            &ListOfExternalStateVariables::addExternalStateVariableOnGroupOfCells<
                ConcreteDryingExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnMesh",
            &ListOfExternalStateVariables::addExternalStateVariableOnMesh<
                TotalFluidPressureExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnGroupOfCells",
            &ListOfExternalStateVariables::addExternalStateVariableOnGroupOfCells<
                TotalFluidPressureExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnMesh",
            &ListOfExternalStateVariables::addExternalStateVariableOnMesh<
                VolumetricDeformationExternalStateVariablePtr > );
    c3.def( "addExternalStateVariableOnGroupOfCells",
            &ListOfExternalStateVariables::addExternalStateVariableOnGroupOfCells<
                VolumetricDeformationExternalStateVariablePtr > );
};
