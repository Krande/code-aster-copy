/**
 * @file ListOfLoadsInterface.cxx
 * @brief Interface python de ListOfLoads
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

#include "PythonBindings/ListOfLoadsInterface.h"
#include "PythonBindings/LoadUtilities.h"
#include "PythonBindings/factory.h"
#include <boost/python.hpp>

namespace py = boost::python;

void exportListOfLoadsToPython() {

    py::class_< ListOfLoads, ListOfLoadsPtr, py::bases< DataStructure > > c1( "ListOfLoads",
                                                                              py::no_init );
    c1.def( "__init__", py::make_constructor( &initFactoryPtr< ListOfLoads > ) );
    c1.def( "__init__", py::make_constructor( &initFactoryPtr< ListOfLoads, ModelPtr > ) );
    c1.def( "isEmpty", &ListOfLoads::isEmpty, R"(
            The list of loads is empty or not.

            Return:
                bool : True if empty
        )",
            ( py::arg( "self" ) ) );
    c1.def( "hasDirichletBC", &ListOfLoads::hasDirichletBC, R"(
            Dirichlet BCs have been added or not ?

            Return:
                bool : True if Dirichlet BCs have been added
        )",
            ( py::arg( "self" ) ) );
    c1.def( "hasExternalLoad", &ListOfLoads::hasExternalLoad, R"(
            External load (= not Dirichlet BCs) have been added or not ?

            Return:
                bool : True if External load have been added
        )",
            ( py::arg( "self" ) ) );
    c1.def( "getDirichletBCs", &ListOfLoads::getDirichletBCs,
            py::return_value_policy< py::copy_const_reference >(), R"(
Return list of DirichletBCs

Returns:
    ListDiriBC: a list of DirichletBC
        )",
            ( py::arg( "self" ) ) );
    c1.def( "getMechanicalLoadsReal", &ListOfLoads::getMechanicalLoadsReal,
            py::return_value_policy< py::copy_const_reference >(), R"(
Return list of real mechanical loads

Returns:
    ListMecaLoadReal: a list of real mechanical loads
        )",
            ( py::arg( "self" ) ) );
    c1.def( "getMechanicalLoadsFunction", &ListOfLoads::getMechanicalLoadsFunction,
            py::return_value_policy< py::copy_const_reference >(), R"(
Return list of Function mechanical loads

Returns:
    ListMecaLoadFunction: a list of Function mechanical loads
        )",
            ( py::arg( "self" ) ) );
#ifdef ASTER_HAVE_MPI
    c1.def( "getParallelMechanicalLoadsReal",
            &ListOfLoads::getParallelMechanicalLoadsReal,
            py::return_value_policy< py::copy_const_reference >(), R"(
Return list of real parallel mechanical loads

Returns:
    ListParaMecaLoadReal: a list of real parallel mechanical loads
        )",
            ( py::arg( "self" ) ) );
    c1.def( "getParallelMechanicalLoadsFunction",
            &ListOfLoads::getParallelMechanicalLoadsFunction,
            py::return_value_policy< py::copy_const_reference >(), R"(
Return list of function parallel mechanical loads

Returns:
    ListParaMecaLoadFunction: a list of function parallel mechanical loads
        )",
            ( py::arg( "self" ) ) );
#endif /* ASTER_HAVE_MPI */
};
