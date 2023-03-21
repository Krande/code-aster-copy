/**
 * @file InterspectralMatrix.cxx
 * @brief Implementation de InterspectralMatrix
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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

#include "LinearAlgebra/InterspectralMatrix.h"

InterspectralMatrix::InterspectralMatrix( const std::string name )
        : DataStructure( name, 8, "INTERSPECTRE" ),
          _refe( JeveuxVectorChar16( getName() + ".REFE" ) ),
          _disc( JeveuxVectorReal( getName() + ".DISC" ) ),
          _vale( JeveuxCollectionReal( getName() + ".VALE" ) ),
          _numi( JeveuxVectorLong( getName() + ".NUMI" ) ),
          _numj( JeveuxVectorLong( getName() + ".NUMJ" ) ),
          _numeOrdre( JeveuxVectorLong( getName() + ".NUME_ORDRE" ) ),
          _noei( JeveuxVectorChar8( getName() + ".NOEI" ) ),
          _noej( JeveuxVectorChar8( getName() + ".NOEJ" ) ),
          _cmpi( JeveuxVectorChar8( getName() + ".CMPI" ) ),
          _cmpj( JeveuxVectorChar8( getName() + ".CMPJ" ) ){};

InterspectralMatrix::InterspectralMatrix() : InterspectralMatrix( ResultNaming::getNewResultName() ){};

std::vector< std::string > InterspectralMatrix::toString( const std::vector< JeveuxChar8 >& vc ) {
    std::vector< std::string > vs;
    vs.reserve( vc.size() );
    for ( auto c : vc )
        vs.push_back( trim( c.toString() ) );
    return vs;
}

std::vector< ASTERINTEGER > InterspectralMatrix::getNumI() const {
    if (!_numi->exists())
        return {};
    _numi->updateValuePointer();
    return _numi->toVector();
};

std::vector< ASTERINTEGER > InterspectralMatrix::getNumJ() const {
    if (!_numj->exists())
        return {};
    _numj->updateValuePointer();
    return _numj->toVector();
};

std::vector< std::string > InterspectralMatrix::getNoeI() const {
    if (!_noei->exists())
        return {};
    _noei->updateValuePointer();
    return toString( _noei->toVector() );
};

std::vector< std::string > InterspectralMatrix::getNoeJ() const {
    if (!_noej->exists())
        return {};
    _noej->updateValuePointer();
    return toString( _noej->toVector() );
};

std::vector< std::string > InterspectralMatrix::getCmpI() const {
    if (!_cmpi->exists())
        return {};
    _cmpi->updateValuePointer();
    return toString( _cmpi->toVector() );
};

std::vector< std::string > InterspectralMatrix::getCmpJ() const {
    if (!_cmpj->exists())
        return {};
    _cmpj->updateValuePointer();
    return toString( _cmpj->toVector() );
};
