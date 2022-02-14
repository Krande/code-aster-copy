#ifndef LOADINTERFACE_H_
#define LOADINTERFACE_H_

/**
 * @file LoadUtilities.h
 * @brief Fichier entete de la classe LoadUtilities
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

#include "aster_pybind.h"
#include "astercxx.h"

#include "Functions/Formula.h"
#include "Functions/Function.h"
#include "Functions/Function2D.h"
#include "Loads/AcousticLoad.h"
#include "Loads/DirichletBC.h"
#include "Loads/MechanicalLoad.h"
#include "Loads/ParallelMechanicalLoad.h"
#include "Loads/ThermalLoad.h"

template < class cppclass, typename... Args >
void addDirichletBCToInterface( py::class_< cppclass, Args... > pyclass ) {

    void ( cppclass::*c1 )( const DirichletBCPtr & ) = &cppclass::addLoad;
    void ( cppclass::*c2 )( const DirichletBCPtr &currentLoad, const FunctionPtr &func ) =
        &cppclass::addLoad;
    void ( cppclass::*c3 )( const DirichletBCPtr &currentLoad, const FormulaPtr &func ) =
        &cppclass::addLoad;
    void ( cppclass::*c4 )( const DirichletBCPtr &currentLoad, const Function2DPtr &func ) =
        &cppclass::addLoad;

    pyclass.def( "addDirichletBC", c1 );
    pyclass.def( "addDirichletBC", c2 );
    pyclass.def( "addDirichletBC", c3 );
    pyclass.def( "addDirichletBC", c4 );
};

template < class cppclass, typename... Args >
void addMechanicalLoadToInterface( py::class_< cppclass, Args... > pyclass ) {

    void ( cppclass::*c5 )( const MechanicalLoadRealPtr & ) = &cppclass::addLoad;
    void ( cppclass::*c6 )( const MechanicalLoadRealPtr &currentLoad, const FunctionPtr &func ) =
        &cppclass::addLoad;
    void ( cppclass::*c7 )( const MechanicalLoadRealPtr &currentLoad, const FormulaPtr &func ) =
        &cppclass::addLoad;
    void ( cppclass::*c8 )( const MechanicalLoadRealPtr &currentLoad, const Function2DPtr &func ) =
        &cppclass::addLoad;

    void ( cppclass::*c9 )( const MechanicalLoadFunctionPtr & ) = &cppclass::addLoad;
    void ( cppclass::*c10 )( const MechanicalLoadFunctionPtr &currentLoad,
                             const FunctionPtr &func ) = &cppclass::addLoad;
    void ( cppclass::*c11 )( const MechanicalLoadFunctionPtr &currentLoad,
                             const FormulaPtr &func ) = &cppclass::addLoad;
    void ( cppclass::*c12 )( const MechanicalLoadFunctionPtr &currentLoad,
                             const Function2DPtr &func ) = &cppclass::addLoad;

    void ( cppclass::*c13 )( const MechanicalLoadComplexPtr & ) = &cppclass::addLoad;
    void ( cppclass::*c14 )( const MechanicalLoadComplexPtr &currentLoad,
                             const FunctionPtr &func ) = &cppclass::addLoad;
    void ( cppclass::*c15 )( const MechanicalLoadComplexPtr &currentLoad, const FormulaPtr &func ) =
        &cppclass::addLoad;
    void ( cppclass::*c16 )( const MechanicalLoadComplexPtr &currentLoad,
                             const Function2DPtr &func ) = &cppclass::addLoad;

    pyclass.def( "addLoad", c5 );
    pyclass.def( "addLoad", c6 );
    pyclass.def( "addLoad", c7 );
    pyclass.def( "addLoad", c8 );
    pyclass.def( "addLoad", c9 );
    pyclass.def( "addLoad", c10 );
    pyclass.def( "addLoad", c11 );
    pyclass.def( "addLoad", c12 );
    pyclass.def( "addLoad", c13 );
    pyclass.def( "addLoad", c14 );
    pyclass.def( "addLoad", c15 );
    pyclass.def( "addLoad", c16 );
};

#ifdef ASTER_HAVE_MPI
template < class cppclass, typename... Args >
void addParallelMechanicalLoadToInterface( py::class_< cppclass, Args... > pyclass ) {

    void ( cppclass::*c9 )( const ParallelMechanicalLoadRealPtr & ) = &cppclass::addLoad;
    void ( cppclass::*c10 )( const ParallelMechanicalLoadRealPtr &currentLoad,
                             const FunctionPtr &func ) = &cppclass::addLoad;
    void ( cppclass::*c11 )( const ParallelMechanicalLoadRealPtr &currentLoad,
                             const FormulaPtr &func ) = &cppclass::addLoad;
    void ( cppclass::*c12 )( const ParallelMechanicalLoadRealPtr &currentLoad,
                             const Function2DPtr &func ) = &cppclass::addLoad;

    void ( cppclass::*c13 )( const ParallelMechanicalLoadFunctionPtr & ) = &cppclass::addLoad;
    void ( cppclass::*c14 )( const ParallelMechanicalLoadFunctionPtr &currentLoad,
                             const FunctionPtr &func ) = &cppclass::addLoad;
    void ( cppclass::*c15 )( const ParallelMechanicalLoadFunctionPtr &currentLoad,
                             const FormulaPtr &func ) = &cppclass::addLoad;
    void ( cppclass::*c16 )( const ParallelMechanicalLoadFunctionPtr &currentLoad,
                             const Function2DPtr &func ) = &cppclass::addLoad;

    pyclass.def( "addLoad", c9 );
    pyclass.def( "addLoad", c10 );
    pyclass.def( "addLoad", c11 );
    pyclass.def( "addLoad", c12 );
    pyclass.def( "addLoad", c13 );
    pyclass.def( "addLoad", c14 );
    pyclass.def( "addLoad", c15 );
    pyclass.def( "addLoad", c16 );
};
#endif /* ASTER_HAVE_MPI */

template < class cppclass, typename... Args >
void addThermalLoadToInterface( py::class_< cppclass, Args... > pyclass ) {

    void ( cppclass::*c13 )( const ThermalLoadRealPtr & ) = &cppclass::addLoad;
    void ( cppclass::*c14 )( const ThermalLoadRealPtr &currentLoad, const FunctionPtr &func ) =
        &cppclass::addLoad;
    void ( cppclass::*c15 )( const ThermalLoadRealPtr &currentLoad, const FormulaPtr &func ) =
        &cppclass::addLoad;
    void ( cppclass::*c16 )( const ThermalLoadRealPtr &currentLoad, const Function2DPtr &func ) =
        &cppclass::addLoad;

    void ( cppclass::*c17 )( const ThermalLoadFunctionPtr & ) = &cppclass::addLoad;
    void ( cppclass::*c18 )( const ThermalLoadFunctionPtr &currentLoad, const FunctionPtr &func ) =
        &cppclass::addLoad;
    void ( cppclass::*c19 )( const ThermalLoadFunctionPtr &currentLoad, const FormulaPtr &func ) =
        &cppclass::addLoad;
    void ( cppclass::*c20 )( const ThermalLoadFunctionPtr &currentLoad,
                             const Function2DPtr &func ) = &cppclass::addLoad;

    pyclass.def( "addLoad", c13 );
    pyclass.def( "addLoad", c14 );
    pyclass.def( "addLoad", c15 );
    pyclass.def( "addLoad", c16 );
    pyclass.def( "addLoad", c17 );
    pyclass.def( "addLoad", c18 );
    pyclass.def( "addLoad", c19 );
    pyclass.def( "addLoad", c20 );
};

template < class cppclass, typename... Args >
void addAcousticLoadToInterface( py::class_< cppclass, Args... > pyclass ) {

    void ( cppclass::*c17 )( const AcousticLoadComplexPtr & ) = &cppclass::addLoad;
    void ( cppclass::*c18 )( const AcousticLoadComplexPtr &currentLoad, const FunctionPtr &func ) =
        &cppclass::addLoad;
    void ( cppclass::*c19 )( const AcousticLoadComplexPtr &currentLoad, const FormulaPtr &func ) =
        &cppclass::addLoad;
    void ( cppclass::*c20 )( const AcousticLoadComplexPtr &currentLoad,
                             const Function2DPtr &func ) = &cppclass::addLoad;

    pyclass.def( "addLoad", c17 );
    pyclass.def( "addLoad", c18 );
    pyclass.def( "addLoad", c19 );
    pyclass.def( "addLoad", c20 );
};

#endif /* LOADINTERFACE_H_ */
