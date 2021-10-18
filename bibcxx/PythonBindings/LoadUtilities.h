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
#include "Functions/Formula.h"
#include "Functions/Function.h"
#include "Functions/Function2D.h"
#include "Loads/AcousticLoad.h"
#include "Loads/DirichletBC.h"
#include "Loads/MechanicalLoad.h"
#include "Loads/ParallelMechanicalLoad.h"
#include "Loads/ThermalLoad.h"

template < class first, typename... Args >
void addDirichletBCToInterface( py::class_< first, Args... > myInstance ) {
    typedef first my;

    void ( my::*c1 )( const DirichletBCPtr & ) = &my::addLoad;
    void ( my::*c2 )( const DirichletBCPtr &currentLoad, const FunctionPtr &func ) =
        &my::addLoad;
    void ( my::*c3 )( const DirichletBCPtr &currentLoad, const FormulaPtr &func ) =
        &my::addLoad;
    void ( my::*c4 )( const DirichletBCPtr &currentLoad, const Function2DPtr &func ) =
        &my::addLoad;

    myInstance.def( "addDirichletBC", c1 );
    myInstance.def( "addDirichletBC", c2 );
    myInstance.def( "addDirichletBC", c3 );
    myInstance.def( "addDirichletBC", c4 );
};

template < class first, typename... Args >
void addMechanicalLoadToInterface( py::class_< first, Args... > myInstance ) {
    typedef first my;

    void ( my::*c5 )( const MechanicalLoadRealPtr & ) = &my::addLoad;
    void ( my::*c6 )( const MechanicalLoadRealPtr &currentLoad, const FunctionPtr &func ) =
        &my::addLoad;
    void ( my::*c7 )( const MechanicalLoadRealPtr &currentLoad, const FormulaPtr &func ) =
        &my::addLoad;
    void ( my::*c8 )( const MechanicalLoadRealPtr &currentLoad, const Function2DPtr &func )
        = &my::addLoad;

    void ( my::*c9 )( const MechanicalLoadFunctionPtr & ) = &my::addLoad;
    void ( my::*c10 )( const MechanicalLoadFunctionPtr &currentLoad, const FunctionPtr &func )=
        &my::addLoad;
    void ( my::*c11 )( const MechanicalLoadFunctionPtr &currentLoad, const FormulaPtr &func ) =
        &my::addLoad;
    void ( my::*c12 )( const MechanicalLoadFunctionPtr &currentLoad, const Function2DPtr &func)
        = &my::addLoad;

    myInstance.def( "addLoad", c5 );
    myInstance.def( "addLoad", c6 );
    myInstance.def( "addLoad", c7 );
    myInstance.def( "addLoad", c8 );
    myInstance.def( "addLoad", c9 );
    myInstance.def( "addLoad", c10 );
    myInstance.def( "addLoad", c11 );
    myInstance.def( "addLoad", c12 );
};

#ifdef ASTER_HAVE_MPI
template < class first, typename... Args >
void
addParallelMechanicalLoadToInterface( py::class_< first, Args... > myInstance ) {
    typedef first my;

    void ( my::*c9 )( const ParallelMechanicalLoadRealPtr & ) = &my::addLoad;
    void ( my::*c10 )( const ParallelMechanicalLoadRealPtr &currentLoad,
                            const FunctionPtr &func )=
        &my::addLoad;
    void ( my::*c11 )( const ParallelMechanicalLoadRealPtr &currentLoad,
                            const FormulaPtr &func ) =
        &my::addLoad;
    void ( my::*c12 )( const ParallelMechanicalLoadRealPtr &currentLoad,
                            const Function2DPtr &func)
        = &my::addLoad;

    void ( my::*c13 )( const ParallelMechanicalLoadFunctionPtr & ) = &my::addLoad;
    void ( my::*c14 )( const ParallelMechanicalLoadFunctionPtr &currentLoad,
                            const FunctionPtr &func )=
        &my::addLoad;
    void ( my::*c15 )( const ParallelMechanicalLoadFunctionPtr &currentLoad,
                            const FormulaPtr &func ) =
        &my::addLoad;
    void ( my::*c16 )( const ParallelMechanicalLoadFunctionPtr &currentLoad,
                            const Function2DPtr &func)
        = &my::addLoad;


    myInstance.def( "addLoad", c9 );
    myInstance.def( "addLoad", c10 );
    myInstance.def( "addLoad", c11 );
    myInstance.def( "addLoad", c12 );
    myInstance.def( "addLoad", c13 );
    myInstance.def( "addLoad", c14 );
    myInstance.def( "addLoad", c15 );
    myInstance.def( "addLoad", c16 );
};
#endif /* ASTER_HAVE_MPI */

template < class first, typename... Args >
void addThermalLoadToInterface( py::class_< first, Args... > myInstance ) {
    typedef first my;

    void ( my::*c13 )( const ThermalLoadRealPtr & ) = &my::addLoad;
    void ( my::*c14 )( const ThermalLoadRealPtr &currentLoad, const FunctionPtr &func ) =
        &my::addLoad;
    void ( my::*c15 )( const ThermalLoadRealPtr &currentLoad, const FormulaPtr &func ) =
        &my::addLoad;
    void ( my::*c16 )( const ThermalLoadRealPtr &currentLoad, const Function2DPtr &func )
        = &my::addLoad;


    void ( my::*c17 )( const ThermalLoadFunctionPtr & ) = &my::addLoad;
    void ( my::*c18 )( const ThermalLoadFunctionPtr &currentLoad, const FunctionPtr &func ) =
        &my::addLoad;
    void ( my::*c19 )( const ThermalLoadFunctionPtr &currentLoad, const FormulaPtr &func ) =
        &my::addLoad;
    void ( my::*c20 )( const ThermalLoadFunctionPtr &currentLoad, const Function2DPtr &func )
        = &my::addLoad;

    myInstance.def( "addLoad", c13 );
    myInstance.def( "addLoad", c14 );
    myInstance.def( "addLoad", c15 );
    myInstance.def( "addLoad", c16 );
    myInstance.def( "addLoad", c17 );
    myInstance.def( "addLoad", c18 );
    myInstance.def( "addLoad", c19 );
    myInstance.def( "addLoad", c20 );
};

template < class first, typename... Args >
void addAcousticLoadToInterface( py::class_< first, Args... > myInstance ) {
    typedef first my;

    void ( my::*c17 )( const AcousticLoadComplexPtr & ) = &my::addLoad;
    void ( my::*c18 )( const AcousticLoadComplexPtr &currentLoad, const FunctionPtr &func ) =
        &my::addLoad;
    void ( my::*c19 )( const AcousticLoadComplexPtr &currentLoad, const FormulaPtr &func ) =
        &my::addLoad;
    void ( my::*c20 )( const AcousticLoadComplexPtr &currentLoad, const Function2DPtr &func )
        = &my::addLoad;

    myInstance.def( "addLoad", c17 );
    myInstance.def( "addLoad", c18 );
    myInstance.def( "addLoad", c19 );
    myInstance.def( "addLoad", c20 );
};

#endif /* LOADINTERFACE_H_ */
