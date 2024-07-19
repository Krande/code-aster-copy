/**
 * @file FieldOnNodes.cxx
 * @brief Implementation de FieldOnNodes vide car FieldOnNodes est un template
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

#include "DataFields/FieldOnNodes.h"

template <>
FieldOnNodesReal FieldOnNodesReal::transform( py::object &func ) const {
    if ( !PyCallable_Check( func.ptr() ) )
        raiseAsterError( "Input parameter to the transform "
                         "method should be a callable Python object" );

    FieldOnNodesReal tmp( *this );
    updateValuePointers();

    ASTERINTEGER size = _values->size();
    for ( auto i = 0; i < size; i++ ) {
        PyObject *res = PyObject_CallFunction( func.ptr(), "d", ( *_values )[i] );
        if ( PyFloat_Check( res ) ) {
            tmp[i] = (ASTERDOUBLE)PyFloat_AsDouble( res );
        } else if ( PyLong_Check( res ) ) {
            tmp[i] = (ASTERDOUBLE)PyLong_AsDouble( res );
        } else {
            raiseAsterError( "Invalid function return type. Expected ASTERDOUBLE." );
        }
        Py_XDECREF( res );
    }
    return tmp;
};

template <>
FieldOnNodesComplex FieldOnNodesComplex::transform( py::object &func ) const {
    if ( !PyCallable_Check( func.ptr() ) )
        raiseAsterError( "Input parameter to the transform "
                         "method should be a callable Python object" );

    FieldOnNodesComplex tmp( *this );
    _values->updateValuePointer();

    ASTERINTEGER size = _values->size();

    Py_complex val;
    for ( auto i = 0; i < size; i++ ) {
        val.real = ( *_values )[i].real();
        val.imag = ( *_values )[i].imag();
        PyObject *res = PyObject_CallFunction( func.ptr(), "D", val );
        if ( PyComplex_Check( res ) ) {
            ASTERDOUBLE re = (ASTERDOUBLE)PyComplex_RealAsDouble( res );
            ASTERDOUBLE im = (ASTERDOUBLE)PyComplex_ImagAsDouble( res );
            tmp[i] = { re, im };
        } else {
            raiseAsterError( "Invalid function return type. Expected ASTERCOMPLEX." );
        }
        Py_XDECREF( res );
    }
    return tmp;
};
