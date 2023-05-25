/**
 * @file FieldOnCells.cxx
 * @brief Implementation de FieldOnCells vide car FieldOnCells est un template
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

#include "DataFields/FieldOnCells.h"

#include "Discretization/ElementaryCharacteristics.h"

// explicit declaration
template <>
FieldOnCells< ASTERDOUBLE >::FieldOnCells( const FiniteElementDescriptorPtr FEDesc,
                                           const std::string &loc, const std::string &quantity,
                                           const BehaviourPropertyPtr behaviour,
                                           const ElementaryCharacteristicsPtr carael )
    : FieldOnCells< ASTERDOUBLE >() {
    std::string inName = getName();
    std::string option;
    std::string nompar;

    if ( loc == "ELGA" ) {
        option = "TOU_INI_ELGA";
    } else if ( loc == "ELNO" ) {
        option = "TOU_INI_ELNO";
    } else {
        AS_ASSERT( false )
    };

    if ( quantity == "SIEF_R" ) {
        nompar = "PSIEF_R";
    } else if ( quantity == "VARI_R" ) {
        nompar = "PVARI_R";
    } else {
        AS_ASSERT( false )
    };

    ASTERINTEGER iret = 0;

    setDescription( FEDesc );

    std::string dcel = " ";

    if ( behaviour ) {
        auto compor = behaviour->getBehaviourField();
        const auto comporName = compor->getName();

        std::string carele = " ";

        if ( carael )
            carele = carael->getName();

        _DCEL = std::make_shared< SimpleFieldOnCellsLong >( inName );
        CALLO_CESVAR( carele, comporName, _dofDescription->getName(), _DCEL->getName() );
        dcel = _DCEL->getName();
    }

    CALLO_ALCHML( _dofDescription->getName(), option, nompar, JeveuxMemoryTypesNames[Permanent],
                  getName(), &iret, dcel );
    AS_ASSERT( iret == 0 );

    updateValuePointers();
};
