/**
 * @file ListOfLoadsInterface.cxx
 * @brief Interface python de ListOfLoads
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

#include "PythonBindings/ListOfLoadsInterface.h"

#include "aster_pybind.h"

#include "PythonBindings/FieldOnNodesInterface.h"
#include "PythonBindings/LoadUtilities.h"

void exportListOfLoadsToPython( py::module_ &mod ) {

    py::class_< ListOfLoads, ListOfLoadsPtr, DataStructure > c1( mod, "ListOfLoads" );
    c1.def( py::init( &initFactoryPtr< ListOfLoads > ) );
    c1.def( py::init( &initFactoryPtr< ListOfLoads, std::string > ) );
    c1.def( py::init( &initFactoryPtr< ListOfLoads, ModelPtr > ) );
    c1.def( py::init( &initFactoryPtr< ListOfLoads, std::string, ModelPtr > ) );
    c1.def( "isBuilt", &ListOfLoads::isBuilt, R"(
            The list of loads has been built or not.

            Returns:
                bool: True if has been built already.
        )" );
    c1.def( "hasDirichletBC", &ListOfLoads::hasDirichletBC, R"(
            Dirichlet BCs have been added or not ?

            Returns:
                bool: True if Dirichlet BCs have been added
        )" );
    c1.def( "hasExternalLoad", &ListOfLoads::hasExternalLoad, R"(
            External load (= not Dirichlet BCs) have been added or not ?

            Returns:
                bool: True if External load have been added
        )" );
    c1.def( "getDirichletBCs", &ListOfLoads::getDirichletBCs, R"(
Return list of DirichletBCs

Returns:
    ListDiriBC: a list of DirichletBC
        )" );
    c1.def( "getMechanicalLoadsReal", &ListOfLoads::getMechanicalLoadsReal,
            R"(
Return list of real mechanical loads

Returns:
    ListMecaLoadReal: a list of real mechanical loads
        )" );
    c1.def( "getMechanicalLoadsFunction", &ListOfLoads::getMechanicalLoadsFunction,
            R"(
Return list of Function mechanical loads

Returns:
    ListMecaLoadFunction: a list of Function mechanical loads
        )" );
#ifdef ASTER_HAVE_MPI
    c1.def( "getParallelMechanicalLoadsReal", &ListOfLoads::getParallelMechanicalLoadsReal,
            R"(
Return list of real parallel mechanical loads

Returns:
    ListParaMecaLoadReal: a list of real parallel mechanical loads
        )" );
    c1.def( "getParallelMechanicalLoadsFunction", &ListOfLoads::getParallelMechanicalLoadsFunction,
            R"(
Return list of function parallel mechanical loads

Returns:
    ListParaMecaLoadFunction: a list of function parallel mechanical loads
        )" );
#endif /* ASTER_HAVE_MPI */

    c1.def( "getModel", &ListOfLoads::getModel, R"(
Return the model used

Returns:
    Model: model used
        )" );
    c1.def( "setDifferentialDisplacement", &ListOfLoads::setDifferentialDisplacement, R"(
Set differential displacement field for DIDI loads
    )" );
    c1.def( "hasDifferential", &ListOfLoads::hasDifferential, R"(
Return True if there are DIDI loads or DIDI Dirichlet BCs
    )" );
    c1.def( "getLoadNames", &ListOfLoads::getLoadNames, R"(
Returns list of load's names.

Returns:
    list[str]: list of load's names.
            )" );
    addDirichletBCToInterface( c1 );
    addMechanicalLoadToInterface( c1 );
    addThermalLoadToInterface( c1 );
    addAcousticLoadToInterface( c1 );
#ifdef ASTER_HAVE_MPI
    addParallelMechanicalLoadToInterface( c1 );
    addParallelThermalLoadToInterface( c1 );
#endif
};
