/**
 * @file CommandSyntax.cxx
 * @brief Implementation of API to CommandSyntax Python object.
 * @section LICENCE
 * Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

#include "Supervis/CommandSyntax.h"

#include "Python.h"
#include "aster.h"
#include "shared_vars.h"

#include "Utilities/CapyConvertibleValue.h"
#include "Utilities/SyntaxDictionary.h"

PyObject *CommandSyntax::py = NULL;

void _check_py() {
    if ( CommandSyntax::py == NULL ) {
        CommandSyntax::py = GetJdcAttr( (char *)"syntax" );
    }
}

CommandSyntax::CommandSyntax( const std::string name ) : _commandName( name ) {
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

void CommandSyntax::define( SyntaxMapContainer &syntax ) {
    PyObject *keywords = syntax.convertToPythonDictionnary();
    PyObject *res = PyObject_CallMethod( _pySyntax, (char *)"define", (char *)"O", keywords );
    if ( res == NULL ) {
        throw std::runtime_error( "Error calling `CommandSyntax.define`." );
    }
    Py_DECREF( res );
    Py_DECREF( keywords );
}

void CommandSyntax::define( py::object keywords ) {
    PyObject *res = PyObject_CallMethod( _pySyntax, (char *)"define", (char *)"O", keywords.ptr() );
    if ( res == NULL ) {
        throw std::runtime_error( "Error calling `CommandSyntax.define`." );
    }
    Py_DECREF( res );
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
