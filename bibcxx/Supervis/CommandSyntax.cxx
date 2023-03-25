/**
 * @file CommandSyntax.cxx
 * @brief Implementation of API to CommandSyntax Python object.
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

#include "aster.h"

#include "Supervis/CommandSyntax.h"

#include "shared_vars.h"

#include "Utilities/CapyConvertibleValue.h"
#include "Utilities/SyntaxDictionary.h"

PyObject *CommandSyntax::py = NULL;

void _check_py() {
    if ( CommandSyntax::py == NULL ) {
        CommandSyntax::py = GetJdcAttr( (char *)"syntax" );
    }
}

CommandSyntax::CommandSyntax( const std::string name ) : _cataPath( name ) {
    _check_py();

    std::string format( "s" );
    _pySyntax = PyObject_CallFunction( CommandSyntax::py, format.c_str(), name.c_str() );
    if ( _pySyntax == NULL ) {
        throw std::runtime_error( "Error during `CommandSyntax.__init__`." );
    }

    register_sh_etape( append_etape( _pySyntax ) );
}

CommandSyntax::~CommandSyntax() { free(); }

void CommandSyntax::free() {
    // already freed?
    if ( _pySyntax == NULL ) {
        return;
    }

    register_sh_etape( pop_etape() );
    PyObject *res = PyObject_CallMethod( _pySyntax, (char *)"free", NULL );
    if ( res == NULL ) {
        throw std::runtime_error( "Error calling `CommandSyntax.free`." );
    }
    Py_DECREF( res );
    Py_CLEAR( _pySyntax );
    Py_CLEAR( py );
}

void CommandSyntax::debugPrint() const {
    PyObject *res = PyObject_CallMethod( _pySyntax, (char *)"__repr__", NULL );
    if ( res == NULL ) {
        throw std::runtime_error( "Error calling `CommandSyntax.__repr__`." );
    }
    Py_DECREF( res );
}

// Syntax checking must be disabled using CapyConvertibleSyntax/SyntaxMapContainer objects,
// because DataStructures are passed as string.
// Using Python keywords arguments, syntax checking should be enabled (defauls to 'true').
void CommandSyntax::define( py::object keywords, bool check_syntax ) {
    // TODO: store _pySyntax as a py::object
    const py::object syntax = py::reinterpret_borrow< py::object >( _pySyntax );
    py::object definePy = syntax.attr( "define" );
    try {
        py::object res = definePy( keywords, py::arg( "check_syntax" ) = check_syntax );
    } catch ( py::error_already_set & ) {
#ifdef ASTER_DEBUG_CXX
        std::abort();
#else
        throw;
#endif
    }
}

void CommandSyntax::define( SyntaxMapContainer &syntax ) {
    const py::object keywords =
        py::reinterpret_steal< py::object >( syntax.convertToPythonDictionnary() );
    define( keywords, false );
}

void CommandSyntax::define( const CapyConvertibleSyntax &syntax ) {
    SyntaxMapContainer container = syntax.toSyntaxMapContainer();
    define( container );
}

void CommandSyntax::setResult( const std::string resultName, const std::string typeSd ) const {
    PyObject *res = PyObject_CallMethod( _pySyntax, (char *)"setResult", (char *)"ss",
                                         resultName.c_str(), typeSd.c_str() );
    if ( res == NULL ) {
        throw std::runtime_error( "Error calling `CommandSyntax.setResult`." );
    }
    Py_DECREF( res );
}
