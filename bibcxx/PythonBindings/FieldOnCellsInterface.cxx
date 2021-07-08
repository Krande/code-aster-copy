/**
 * @file FieldOnCellsInterface.cxx
 * @brief Interface python de FieldOnCells
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include <boost/python.hpp>

namespace py = boost::python;
#include <PythonBindings/factory.h>

#include "PythonBindings/DataStructureInterface.h"
#include "PythonBindings/FieldOnCellsInterface.h"

void exportFieldOnCellsToPython() {
    py::class_<  FieldOnCellsReal, FieldOnCellsRealPtr,
            py::bases< DataField > >( "FieldOnCellsReal", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< FieldOnCellsReal >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr<  FieldOnCellsReal, std::string >) )
        .def(py::init<const FieldOnCellsReal&>() )
        .def( "exportToSimpleFieldOnCells",
              &FieldOnCellsReal::exportToSimpleFieldOnCells )
        .def( "getModel", &FieldOnCellsReal::getModel )
        .def( "setDescription", &FieldOnCellsReal::setDescription )
        .def( "setModel", &FieldOnCellsReal::setModel )
        .def( "build", &FieldOnCellsReal::build )
        .def( "transform", &FieldOnCellsReal::transform<ASTERDOUBLE> )
        .def( "getValues", &FieldOnCellsReal::getValues,
              py::return_value_policy<py::copy_const_reference>())
        .def( "__getitem__",
              +[]( const FieldOnCellsReal& v, ASTERINTEGER i ) { return v[i]; } )
        .def( "__setitem__",
              +[]( FieldOnCellsReal &v, ASTERINTEGER i, float f ) { return v.operator[]( i )=f; } )
        .def( "__len__",
              +[]( const FieldOnCellsReal& v ) { return v.size(); } )
        .def( py::self + py::self)
        .def( py::self - py::self)
        .def( py::self += py::self )
        .def( py::self -= py::self )
        .def( py::self * float() )
        .def( float() * py::self )
        .def( - py::self )
        .def( "printMedFile", &FieldOnCellsReal::printMedFile, R"(
                  Print the field in MED format.

                  Arguments:
                  filename (str): Path to the file to be printed.

                  Returns:
                  bool: *True* if succeeds, *False* otherwise.
                        )",
                              ( py::arg( "self" ), py::arg( "filename" ) )  )
        .def( "norm", &FieldOnCellsReal::norm<ASTERDOUBLE>,
               R"(
                  Return the euclidean norm of the field

                  Argument:
                  normType: "NORM_1", "NORM_2", "NORM_INFINITY"

                  Returns:
                  double: euclidean norm
                        )" );


    
    py::class_< FieldOnCellsComplex, FieldOnCellsComplexPtr,
            py::bases< DataField > >( "FieldOnCellsComplex", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< FieldOnCellsComplex >) )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< FieldOnCellsComplex, std::string >) )
        .def(py::init<const FieldOnCellsComplex&>() )
        .def( "getModel", &FieldOnCellsComplex::getModel )
        .def( "setDescription", &FieldOnCellsComplex::setDescription )
        .def( "setModel", &FieldOnCellsComplex::setModel )
        .def( "build", &FieldOnCellsComplex::build )
        .def( "transform", &FieldOnCellsComplex::transform<ASTERCOMPLEX> )
        .def( "getValues", &FieldOnCellsReal::getValues,
               py::return_value_policy<py::copy_const_reference>())
        .def( "__getitem__",
              +[]( const FieldOnCellsComplex& v, int i ) { return v[i]; } )
        .def( "__setitem__",
              +[]( FieldOnCellsComplex &v, ASTERINTEGER i, ASTERCOMPLEX f )
               { return v.operator[]( i )=f; } )
        .def( "__len__",
              +[]( const FieldOnCellsReal& v ) { return v.size(); } )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( py::self += py::self )
        .def( py::self -= py::self )
        .def( py::self * float() )
        .def( float() * py::self )
        .def( "printMedFile", &FieldOnCellsComplex::printMedFile, R"(
                  Print the field in MED format.

                  Arguments:
                  filename (str): Path to the file to be printed.

                  Returns:
                  bool: *True* if succeeds, *False* otherwise.
                        )",
                  ( py::arg( "self" ), py::arg( "filename" ) )  );
};
