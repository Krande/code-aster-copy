/**
 * @file ResultNaming.cxx
 * @brief Implementation of automatic naming of jeveux objects.
 * @section LICENCE
 * Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
 * This file is part of code_aster.
 *
 * code_aster is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * code_aster is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with code_aster.  If not, see <http://www.gnu.org/licenses/>.

 * person_in_charge: mathieu.courtois@edf.fr
 */

#include "astercxx.h"

#include "Functions/Formula.h"

#include "aster_pybind.h"

#include "Supervis/ResultNaming.h"
#include "Utilities/Tools.h"

#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

Formula::Formula( const std::string name )
    : GenericFunction( name, "FORMULE", "INTERPRE" ),
      _variables( JeveuxVectorChar24( getName() + ".NOVA" ) ),
      _pointers( JeveuxVectorLong( getName() + ".ADDR" ) ),
      _expression( "" ),
      _code( py::none() ),
      _context( py::dict() ) {}

Formula::Formula() : Formula::Formula( ResultNaming::getNewResultName() ) {
    propertyAllocate();
    _pointers->allocate( 2 );
}

void Formula::setVariables( const VectorString &names ) {
    const int nbvar = names.size();
    _variables->allocate( nbvar );

    VectorString::const_iterator varIt = names.begin();
    int idx = 0;
    for ( ; varIt != names.end(); ++varIt ) {
        ( *_variables )[idx] = *varIt;
        ++idx;
    }
}

VectorString Formula::getVariables() const {
    _variables->updateValuePointer();
    long nbvars = _variables->size();
    VectorString vars;
    for ( int i = 0; i < nbvars; ++i ) {
        vars.push_back( ( *_variables )[i].rstrip() );
    }
    return vars;
}

void Formula::setExpression( const std::string expression ) {
    std::string name( strip( "code_of_" + getName() ) );
    _expression = expression;
    PyCompilerFlags flags;
    flags.cf_flags = CO_FUTURE_DIVISION;

    PyObject *pcode;
    pcode = Py_CompileStringFlags( _expression.c_str(), name.c_str(), Py_eval_input, &flags );
    _code = py::reinterpret_steal< py::object >( pcode );

    _pointers->updateValuePointer();
    ( *_pointers )[0] = (ASTERINTEGER)_code.ptr();
    if ( _code.ptr() == nullptr ) {
        PyErr_Print();
        throw std::runtime_error( "Invalid syntax in expression." );
    }
}

void Formula::setContext( py::object context ) {
    if ( !PyDict_Check( context.ptr() ) ) {
        throw std::runtime_error( "Formula: 'dict' object is expected." );
    }
    _context = context;
    _pointers->updateValuePointer();
    ( *_pointers )[1] = (ASTERINTEGER)_context.ptr();
}

VectorString Formula::getProperties() const {
    _property->updateValuePointer();
    VectorString prop;
    for ( int i = 0; i < 6; ++i ) {
        prop.push_back( ( *_property )[i].rstrip() );
    }
    return prop;
}

VectorReal Formula::evaluate( const VectorReal &values ) const {
    int iret = 0;
    VectorString vars = getVariables();
    VectorReal result = evaluate_formula( _code, _context, vars, values, &iret );
    if ( iret == 1 ) {
        const long nbvars = vars.size();
        const long nbvalues = values.size();
        throw std::runtime_error( "Expecting exactly " + std::to_string( nbvars ) +
                                  " values, not " + std::to_string( nbvalues ) );
    } else if ( iret == 4 ) {
        throw std::runtime_error( "Formula: Error during evaluation." );
    }
    return result;
}

/* functions shared with evaluation from Fortran */
VectorReal evaluate_formula( const py::object &code, const py::object &globals,
                             const VectorString &variables, const VectorReal &values,
                             int *retcode ) {
    if ( code.is_none() ) {
        std::cerr << "Formula has no expression:" << std::endl;
        *retcode = 4;
        return VectorReal( 0., 0 );
    }
    const long nbvars = variables.size();
    const long nbvalues = values.size();
    if ( nbvalues != nbvars ) {
        *retcode = 1;
        return VectorReal( 0., 0 );
    }

    py::object locals = py::dict();
    for ( int i = 0; i < nbvars; ++i ) {
        locals[strip( variables[i] ).c_str()] = values[i];
    }

    // res = py::eval( expression, globals, locals );
    PyObject *res = PyEval_EvalCode( code.ptr(), globals.ptr(), locals.ptr() );
    if ( res == NULL ) {
        std::cerr << "Evaluation failed with: ";
        PyObject_Print( locals.ptr(), stderr, 0 );
        std::cerr << std::endl;
        if ( PyErr_Occurred() ) {
            std::cerr << "Detailed traceback of evaluation:" << std::endl;
            PyErr_Print();
            throw py::error_already_set();
        }
        *retcode = 4;
        return VectorReal( 0., 0 );
    }

    VectorReal result;
    if ( PyTuple_Check( res ) ) {
        py::tuple tup = py::reinterpret_borrow< py::tuple >( res );
        for ( auto &value : tup ) {
            result.push_back( value.cast< double >() );
        }
    } else if ( PyComplex_Check( res ) ) {
        result.push_back( PyComplex_RealAsDouble( res ) );
        result.push_back( PyComplex_ImagAsDouble( res ) );
    } else {
        result.push_back( PyFloat_AsDouble( res ) );
    }
    Py_DECREF( res );

    return result;
}

/* Interface for Fortran calls */
void DEFPPSPPPPP( EVAL_FORMULA, eval_formula, ASTERINTEGER *pcode, ASTERINTEGER *pglobals,
                  char *array_vars, STRING_SIZE lenvars, ASTERDOUBLE *array_values,
                  ASTERINTEGER *nbvar, _OUT ASTERINTEGER *iret, _IN ASTERINTEGER *nbres,
                  _OUT ASTERDOUBLE *result ) {
    const py::object code = py::reinterpret_borrow< py::object >( (PyObject *)( *pcode ) );
    const py::object globals = py::reinterpret_borrow< py::object >( (PyObject *)( *pglobals ) );

    VectorString vars;
    VectorReal values;
    for ( int i = 0; i < *nbvar; ++i ) {
        vars.push_back( std::string( array_vars + i * lenvars, lenvars ) );
        values.push_back( array_values[i] );
    }

    int ret = 0;
    VectorReal retvalues = evaluate_formula( code, globals, vars, values, &ret );
    *iret = (ASTERINTEGER)ret;
    if ( ret == 0 ) {
        for ( long i = 0; i < int( retvalues.size() ) && i < ( *nbres ); ++i ) {
            result[i] = (ASTERDOUBLE)retvalues[i];
        }
    }
    return;
}
