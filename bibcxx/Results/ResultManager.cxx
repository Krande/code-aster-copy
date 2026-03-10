/**
 * @file ResultManager.cxx
 * @brief Implementation de ResultManager
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2026  EDF www.code-aster.org
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

#include "Results/ResultManager.h"

void ResultManager::addFieldToResult( const std::string &nomSymb, const std::string &name,
                                      ASTERINTEGER storageIndex ) {
    if ( _currentResult != nullptr ) {
        _currentResult->addFieldFromString( nomSymb, name, storageIndex );
    }
};

extern "C" void DEFSSSP( ADD_FIELD_IN_CURRENT_RESULT, add_field_in_current_result, _IN char *nomres,
                         _IN STRING_SIZE lres, _IN char *nomcham, _IN STRING_SIZE lnoch,
                         _IN char *nomsym, _IN STRING_SIZE lnosy, ASTERINTEGER *iordr ) {
    auto instance = ResultManager::getInstance();
    auto resu = ResultManager::getInstance()->getCurrentResult();
    std::string resuName( nomres, lres ), test1( nomcham, lnoch ), test2( nomsym, lnosy );
    auto name = strip( test2 );
    auto nomSymb = strip( test1 );
    instance->addFieldToResult( nomSymb, name, *iordr );
};
