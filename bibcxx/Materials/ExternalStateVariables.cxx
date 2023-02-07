/**
 * @file ExternalStateVariables.cxx
 * @brief Implementation of ExternalStateVariables
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
#include "Materials/ExternalStateVariables.h"

#include "astercxx.h"

void ExternalStateVariable::setReferenceValue( const ASTERDOUBLE &value ) {
    AS_ASSERT( ExternalVariableTraits::externVarHasRefeValue( _type ) );
    _refValue = value;
};

void EvolutionParameter::setLeftExtension( const std::string typeExtension ) {
    if ( typeExtension == "EXCLU" || typeExtension == "CONSTANT" || typeExtension == "LINEAIRE" ) {
        _leftExtension = typeExtension;
    } else {
        AS_ABORT( "Unknown type of extension for function (left)" );
    };
};
void EvolutionParameter::setRightExtension( const std::string typeExtension ) {
    if ( typeExtension == "EXCLU" || typeExtension == "CONSTANT" || typeExtension == "LINEAIRE" ) {
        _rightExtension = typeExtension;
    } else {
        AS_ABORT( "Unknown type of extension for function (right)" );
    };
};

/** @brief Get transient result */
TransientResultPtr ExternalStateVariable::getTransientResult() const {
    if ( _evolParameter ) {
        return _evolParameter->getTransientResult();
    }
    return nullptr;
}
