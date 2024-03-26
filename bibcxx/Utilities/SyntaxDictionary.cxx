/**
 * @file SyntaxDictionary.cxx
 * @brief Implementation de SyntaxMapContainer
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

#include "Utilities/SyntaxDictionary.h"

PyObject *SyntaxMapContainer::convertToPythonDictionnary( PyObject *returnDict ) const {
    if ( returnDict == NULL )
        returnDict = PyDict_New();

    for ( auto &[key, var] : container ) {
        if ( std::holds_alternative< ASTERINTEGER >( var ) ) {
            const ASTERINTEGER &tmp = std::get< ASTERINTEGER >( var );
            PyObject *value = PyLong_FromLong( tmp );
            PyDict_SetItemString( returnDict, key.c_str(), value );
            Py_DECREF( value );
        } else if ( std::holds_alternative< VectorLong >( var ) ) {
            const VectorLong &currentList = std::get< VectorLong >( var );
            PyObject *listValues = PyList_New( currentList.size() );
            int count = 0;
            for ( auto &val : currentList ) {
                PyList_SetItem( listValues, count, PyLong_FromLong( val ) );
                ++count;
            }
            PyDict_SetItemString( returnDict, key.c_str(), listValues );
            Py_DECREF( listValues );
        } else if ( std::holds_alternative< std::string >( var ) ) {
            const std::string &tmp = std::get< std::string >( var );
            PyObject *value = PyUnicode_FromString( tmp.c_str() );
            PyDict_SetItemString( returnDict, key.c_str(), value );
            Py_DECREF( value );
        } else if ( std::holds_alternative< VectorString >( var ) ) {
            const VectorString &currentList = std::get< VectorString >( var );
            PyObject *listValues = PyList_New( currentList.size() );
            int count = 0;
            for ( auto &val : currentList ) {
                PyList_SetItem( listValues, count, PyUnicode_FromString( val.c_str() ) );
                ++count;
            }
            PyDict_SetItemString( returnDict, key.c_str(), listValues );
            Py_DECREF( listValues );
        } else if ( std::holds_alternative< ASTERDOUBLE >( var ) ) {
            const ASTERDOUBLE &tmp = std::get< ASTERDOUBLE >( var );
            PyObject *value = PyFloat_FromDouble( tmp );
            PyDict_SetItemString( returnDict, key.c_str(), value );
            Py_DECREF( value );
        } else if ( std::holds_alternative< VectorReal >( var ) ) {
            const VectorReal &currentList = std::get< VectorReal >( var );
            PyObject *listValues = PyList_New( currentList.size() );
            int count = 0;
            for ( auto &val : currentList ) {
                PyList_SetItem( listValues, count, PyFloat_FromDouble( val ) );
                ++count;
            }
            PyDict_SetItemString( returnDict, key.c_str(), listValues );
            Py_DECREF( listValues );
        } else if ( std::holds_alternative< ASTERCOMPLEX >( var ) ) {
            const ASTERCOMPLEX &tmp = std::get< ASTERCOMPLEX >( var );
            PyObject *value = PyComplex_FromDoubles( tmp.real(), tmp.imag() );
            PyDict_SetItemString( returnDict, key.c_str(), value );
            Py_DECREF( value );
        } else if ( std::holds_alternative< VectorComplex >( var ) ) {
            const VectorComplex &currentList = std::get< VectorComplex >( var );
            PyObject *listValues = PyList_New( currentList.size() );
            int count = 0;
            for ( auto &val : currentList ) {
                PyList_SetItem( listValues, count,
                                PyComplex_FromDoubles( val.real(), val.imag() ) );
                ++count;
            }
            PyDict_SetItemString( returnDict, key.c_str(), listValues );
            Py_DECREF( listValues );
        } else if ( std::holds_alternative< ListSyntaxMapContainer >( var ) ) {
            const ListSyntaxMapContainer &tmp = std::get< ListSyntaxMapContainer >( var );
            PyObject *list_F = PyList_New( tmp.size() );
            int count = 0;
            for ( ListSyntaxMapContainerCIter iter2 = tmp.begin(); iter2 != tmp.end(); ++iter2 ) {
                PyObject *currentDict = iter2->convertToPythonDictionnary();
                PyList_SetItem( list_F, count, currentDict );
                ++count;
            }
            PyDict_SetItemString( returnDict, key.c_str(), list_F );
            Py_DECREF( list_F );
        }
    }
    return returnDict;
};

py::dict SyntaxMapContainer::convertToPyDict() const {
    py::dict returnDict;

    for ( auto &[key, var] : container ) {
        if ( std::holds_alternative< ASTERINTEGER >( var ) ) {
            const ASTERINTEGER &tmp = std::get< ASTERINTEGER >( var );
            returnDict[key.c_str()] = tmp;
        } else if ( std::holds_alternative< VectorLong >( var ) ) {
            const VectorLong &currentList = std::get< VectorLong >( var );
            py::list listValues;
            for ( auto &val : currentList ) {
                listValues.append( val );
            }
            returnDict[key.c_str()] = listValues;
        } else if ( std::holds_alternative< std::string >( var ) ) {
            const std::string &tmp = std::get< std::string >( var );
            returnDict[key.c_str()] = tmp;
        } else if ( std::holds_alternative< VectorString >( var ) ) {
            const VectorString &currentList = std::get< VectorString >( var );
            py::list listValues;
            for ( auto &val : currentList ) {
                listValues.append( val );
            }
            returnDict[key.c_str()] = listValues;
        } else if ( std::holds_alternative< ASTERDOUBLE >( var ) ) {
            const ASTERDOUBLE &tmp = std::get< ASTERDOUBLE >( var );
            returnDict[key.c_str()] = tmp;
        } else if ( std::holds_alternative< VectorReal >( var ) ) {
            const VectorReal &currentList = std::get< VectorReal >( var );
            py::list listValues;
            for ( auto &val : currentList ) {
                listValues.append( val );
            }
            returnDict[key.c_str()] = listValues;
        } else if ( std::holds_alternative< ASTERCOMPLEX >( var ) ) {
            const ASTERCOMPLEX &tmp = std::get< ASTERCOMPLEX >( var );
            returnDict[key.c_str()] = tmp;
        } else if ( std::holds_alternative< VectorComplex >( var ) ) {
            const VectorComplex &currentList = std::get< VectorComplex >( var );
            py::list listValues;
            for ( auto &val : currentList ) {
                listValues.append( val );
            }
            returnDict[key.c_str()] = listValues;
        } else if ( std::holds_alternative< ListSyntaxMapContainer >( var ) ) {
            const ListSyntaxMapContainer &tmp = std::get< ListSyntaxMapContainer >( var );
            py::list list_F;
            for ( ListSyntaxMapContainerCIter iter2 = tmp.begin(); iter2 != tmp.end(); ++iter2 ) {
                auto currentDict = iter2->convertToPyDict();
                list_F.append( currentDict );
            }
            returnDict[key.c_str()] = list_F;
        }
    }
    return returnDict;
};

SyntaxMapContainer operator+( const SyntaxMapContainer &toAdd1, const SyntaxMapContainer &toAdd2 ) {
    SyntaxMapContainer retour = toAdd1;
    retour += toAdd2;
    return retour;
};
