/**
 * @file FieldOnCells.cxx
 * @brief Implementation de FieldOnCells vide car FieldOnCells est un template
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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

#include "DataFields/FieldOnCells.h"

template <>
FieldOnCellsReal FieldOnCellsReal::transform( py::object &func ) const {
    if ( !PyCallable_Check( func.ptr() ) )
        raiseAsterError( "Input parameter to the transform "
                         "method should be a callable Python object" );

    FieldOnCellsReal tmp( *this );
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
}

template <>
FieldOnCellsComplex FieldOnCellsComplex::transform( py::object &func ) const {
    if ( !PyCallable_Check( func.ptr() ) )
        raiseAsterError( "Input parameter to the transform "
                         "method should be a callable Python object" );

    FieldOnCellsComplex tmp( *this );
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
}

template <>
ASTERDOUBLE FieldOnCellsReal::norm( const std::string normType ) const {
    AS_ASSERT( normType == "NORM_1" || normType == "NORM_2" || normType == "NORM_INFINITY" );

    ASTERDOUBLE norme = 0.0;
    ASTERINTEGER beg = 0, end = 0, nbgrp = 0;

    const int rank = getMPIRank();

    _values->updateValuePointer();
    _descriptor->updateValuePointer();

    JeveuxVectorLong CellsRank = getMesh()->getCellsOwner();
    CellsRank->updateValuePointer();

    auto collec = _dofDescription->getListOfGroupsOfElements();
    nbgrp = ( *_descriptor )[1];

    for ( auto i = 0; i < nbgrp; i++ ) {

        ASTERINTEGER adress = ( *_descriptor )[4 + i];
        if ( ( *_descriptor )[adress + 2] == 0 )
            continue;

        ASTERINTEGER nel = ( *_descriptor )[adress];
        auto liel = ( *collec )[i + 1];
        liel->updateValuePointer();

        if ( normType == "NORM_1" ) {
            for ( auto p = 0; p < nel; p++ ) {

                if ( ( *CellsRank )[( *liel )[p] - 1] != rank )
                    continue;
                beg = ( *_descriptor )[adress + 3 + 4 * p + 4] - 1;
                end = beg + ( *_descriptor )[adress + 3 + 4 * p + 3];

                for ( int pos = beg; pos < end; ++pos ) {
                    norme += std::abs( ( *this )[pos] );
                }
            }
        } else if ( normType == "NORM_2" ) {
            for ( auto p = 0; p < nel; p++ ) {

                if ( ( *CellsRank )[( *liel )[p] - 1] != rank )
                    continue;
                beg = ( *_descriptor )[adress + 3 + 4 * p + 4] - 1;
                end = beg + ( *_descriptor )[adress + 3 + 4 * p + 3];

                for ( int pos = beg; pos < end; ++pos ) {
                    norme += ( *this )[pos] * ( *this )[pos];
                }
            }
        } else if ( normType == "NORM_INFINITY" ) {
            for ( auto p = 0; p < nel; p++ ) {

                if ( ( *CellsRank )[( *liel )[p] - 1] != rank )
                    continue;
                beg = ( *_descriptor )[adress + 3 + 4 * p + 4] - 1;
                end = beg + ( *_descriptor )[adress + 3 + 4 * p + 3];

                for ( int pos = beg; pos < end; ++pos ) {
                    norme = std::max( norme, std::abs( ( *this )[pos] ) );
                }
            }
        }
    }

#ifdef ASTER_HAVE_MPI
    if ( getMesh()->isParallel() ) {
        ASTERDOUBLE norm2 = norme;
        if ( normType == "NORM_1" || normType == "NORM_2" )
            AsterMPI::all_reduce( norm2, norme, MPI_SUM );
        else
            AsterMPI::all_reduce( norm2, norme, MPI_MAX );
    }
#endif

    if ( normType == "NORM_2" )
        norme = std::sqrt( norme );

    return norme;
}

template <>
ASTERDOUBLE FieldOnCellsReal::dot( const FieldOnCellsPtr &tmp ) const {
    tmp->updateValuePointers();
    _values->updateValuePointer();
    ASTERINTEGER taille = _values->size();

    if ( taille != tmp->size() )
        raiseAsterError( "Incompatible size" );

    ASTERDOUBLE ret = 0.0;
    for ( auto pos = 0; pos < taille; ++pos ) {
        ret += ( *this )[pos] * ( *tmp )[pos];
    }
    return ret;
};

template <>
void FieldOnCellsReal::checkInternalStateVariables(
    const ConstantFieldOnCellsChar16Ptr prevBehaviour,
    const ConstantFieldOnCellsChar16Ptr currBehaviour,
    const FiniteElementDescriptorPtr newFEDesc ) {

    std::string prevBehaviourStr = " ";
    std::string currBehaviourStr = " ";
    std::string variStr = " ";
    std::string FEDescStr = " ";

    currBehaviourStr = currBehaviour->getName();
    variStr = this->getName();
    FEDescStr = newFEDesc->getName();

    if ( prevBehaviour != nullptr ) {
        prevBehaviourStr = prevBehaviour->getName();
    };

    CALL_CHCKVARI( prevBehaviourStr, currBehaviourStr, variStr, FEDescStr );
}

template <>
ASTERINTEGER FieldOnCellsReal::compareShape( const FieldOnCellsRealPtr fieldModel,
                                             const bool projectOnLigrel,
                                             const std::string paraName ) {

    ASTERINTEGER iret = 0;
    std::string fieldModelStr = ljust( strip( fieldModel->getName() ), 24, ' ' );
    std::string fieldStr = ljust( strip( this->getName() ), 24, ' ' );
    std::string paraNameStr = ljust( strip( paraName ), 24, ' ' );

    CALL_COMPAREFIELDSHAPE( fieldModelStr, fieldStr, (ASTERLOGICAL *)&projectOnLigrel, paraNameStr,
                            &iret );

    return iret;
}
