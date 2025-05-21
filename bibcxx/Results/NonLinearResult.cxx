/**
 * @file NonLinearResult.cxx
 * @brief Implementation de NonLinearResult
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

#include "Results/NonLinearResult.h"

#include "Behaviour_type.h"

#include "Supervis/Exceptions.h"

#include <filesystem>

JeveuxVectorReal NonLinearResult::_mata( "PYTHON.TANGENT.MATA" );
JeveuxVectorReal NonLinearResult::_matc( "PYTHON.TANGENT.MATC" );

void NonLinearResult::setContact( const ContactPtr contact, const ASTERINTEGER &rank ) {
    if ( !contact )
        raiseAsterError( "ValueError: Contact is empty" );

    _mapContact[rank] = contact;
    const auto fed = contact->getFiniteElementDescriptor();
    _fieldBuilder.addFiniteElementDescriptor( fed );
};

void NonLinearResult::setContact( const ContactPtr contact ) {
    auto indexes = getIndexes();
    for ( auto &index : indexes ) {
        setContact( contact, index );
    }
};

VectorReal NonLinearResult::getTangentMatrix( const std::string &suffix ) {
    if ( suffix == "MATA" && _mata.exists() ) {
        _mata->updateValuePointer();
        return _mata->toVector();
    } else if ( suffix == "MATC" && _matc.exists() ) {
        _matc->updateValuePointer();
        return _matc->toVector();
    } else
        return {};
};

void NonLinearResult::printMedFile( const std::filesystem::path &fileName, std::string medName,
                                    bool local, bool internalVar ) const {
    const auto indexes = getIndexes();
    bool mFront = false;
    if ( _dictOfMapOfConstantFieldOnCellsChar16.find( "COMPORTEMENT" ) !=
         _dictOfMapOfConstantFieldOnCellsChar16.end() ) {
        for ( const auto &index : indexes ) {
            const auto compor = getConstantFieldOnCellsChar16( "COMPORTEMENT", index, true );
            const auto size = compor->size();
            const auto values = compor->getValues();
            for ( const auto &tmp : values ) {
                const auto &val2 = tmp.getValues();
                // MGIS_ADDR-1 because MGIS_ADDR come from Fortran
                const auto &val23 = strip( val2[MGIS_ADDR - 1] );
                if ( val23 != "VIDE" && val23 != "" ) {
                    mFront = true;
                    break;
                }
            }
            if ( mFront )
                break;
            break;
        }
    }
    Result::printMedFile( fileName, medName, local, !mFront );
    return;
}
