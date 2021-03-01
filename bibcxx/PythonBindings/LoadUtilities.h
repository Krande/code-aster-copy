#ifndef LOADINTERFACE_H_
#define LOADINTERFACE_H_

/**
 * @file LoadUtilities.h
 * @brief Fichier entete de la classe LoadUtilities
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

#include "astercxx.h"

#include <boost/python.hpp>

namespace py = boost::python;
#include "Loads/DirichletBC.h"
#include "Loads/MechanicalLoad.h"
#include "Loads/AcousticLoad.h"
#include "Loads/ThermalLoad.h"
#include "Loads/ParallelMechanicalLoad.h"

template < class firstClass, typename... Args >
void addDirichletBCToInterface( py::class_< firstClass, Args... > myInstance ) {
    typedef firstClass myClass;

    void ( myClass::*c1 )( const DirichletBCPtr & ) = &myClass::addLoad;
    void ( myClass::*c2 )( const DirichletBCPtr &currentLoad, const FunctionPtr &func ) =
        &myClass::addLoad;
    void ( myClass::*c3 )( const DirichletBCPtr &currentLoad, const FormulaPtr &func ) =
        &myClass::addLoad;
    void ( myClass::*c4 )( const DirichletBCPtr &currentLoad, const Function2DPtr &func ) =
        &myClass::addLoad;

    myInstance.def( "addDirichletBC", c1 );
    myInstance.def( "addDirichletBC", c2 );
    myInstance.def( "addDirichletBC", c3 );
    myInstance.def( "addDirichletBC", c4 );
};

template < class firstClass, typename... Args >
void addMechanicalLoadToInterface( py::class_< firstClass, Args... > myInstance ) {
    typedef firstClass myClass;

    void ( myClass::*c5 )( const GenericMechanicalLoadPtr & ) = &myClass::addLoad;
    void ( myClass::*c6 )( const GenericMechanicalLoadPtr &currentLoad, const FunctionPtr &func ) =
        &myClass::addLoad;
    void ( myClass::*c7 )( const GenericMechanicalLoadPtr &currentLoad, const FormulaPtr &func ) =
        &myClass::addLoad;
    void ( myClass::*c8 )( const GenericMechanicalLoadPtr &currentLoad, const Function2DPtr &func )
        = &myClass::addLoad;

    myInstance.def( "addLoad", c5 );
    myInstance.def( "addLoad", c6 );
    myInstance.def( "addLoad", c7 );
    myInstance.def( "addLoad", c8 );
};

#ifdef ASTER_HAVE_MPI
template < class firstClass, typename... Args >
void
addParallelMechanicalLoadToInterface( py::class_< firstClass, Args... > myInstance ) {
    typedef firstClass myClass;

    void ( myClass::*c9 )( const ParallelMechanicalLoadPtr & ) = &myClass::addLoad;
    void ( myClass::*c10 )( const ParallelMechanicalLoadPtr &currentLoad, const FunctionPtr &func )=
        &myClass::addLoad;
    void ( myClass::*c11 )( const ParallelMechanicalLoadPtr &currentLoad, const FormulaPtr &func ) =
        &myClass::addLoad;
    void ( myClass::*c12 )( const ParallelMechanicalLoadPtr &currentLoad, const Function2DPtr &func)
        = &myClass::addLoad;

    myInstance.def( "addLoad", c9 );
    myInstance.def( "addLoad", c10 );
    myInstance.def( "addLoad", c11 );
    myInstance.def( "addLoad", c12 );
};
#endif /* ASTER_HAVE_MPI */

template < class firstClass, typename... Args >
void addThermalLoadToInterface( py::class_< firstClass, Args... > myInstance ) {
    typedef firstClass myClass;

    void ( myClass::*c13 )( const ThermalLoadPtr & ) = &myClass::addLoad;
    void ( myClass::*c14 )( const ThermalLoadPtr &currentLoad, const FunctionPtr &func ) =
        &myClass::addLoad;
    void ( myClass::*c15 )( const ThermalLoadPtr &currentLoad, const FormulaPtr &func ) =
        &myClass::addLoad;
    void ( myClass::*c16 )( const ThermalLoadPtr &currentLoad, const Function2DPtr &func )
        = &myClass::addLoad;

    myInstance.def( "addLoad", c13 );
    myInstance.def( "addLoad", c14 );
    myInstance.def( "addLoad", c15 );
    myInstance.def( "addLoad", c16 );
};

template < class firstClass, typename... Args >
void addAcousticLoadToInterface( py::class_< firstClass, Args... > myInstance ) {
    typedef firstClass myClass;

    void ( myClass::*c17 )( const AcousticLoadPtr & ) = &myClass::addLoad;
    void ( myClass::*c18 )( const AcousticLoadPtr &currentLoad, const FunctionPtr &func ) =
        &myClass::addLoad;
    void ( myClass::*c19 )( const AcousticLoadPtr &currentLoad, const FormulaPtr &func ) =
        &myClass::addLoad;
    void ( myClass::*c20 )( const AcousticLoadPtr &currentLoad, const Function2DPtr &func )
        = &myClass::addLoad;

    myInstance.def( "addLoad", c17 );
    myInstance.def( "addLoad", c18 );
    myInstance.def( "addLoad", c19 );
    myInstance.def( "addLoad", c20 );
};

#endif /* LOADINTERFACE_H_ */
