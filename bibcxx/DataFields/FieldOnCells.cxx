/**
 * @file FieldOnCells.cxx
 * @brief Implementation de FieldOnCells vide car FieldOnCells est un template
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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
FieldOnCells< ASTERDOUBLE >::FieldOnCells( const ModelPtr &model,
                                           const BehaviourPropertyPtr behaviour,
                                           const std::string &typcham,
                                           const ElementaryCharacteristicsPtr carael )
    : FieldOnCells< ASTERDOUBLE >() {
    std::string inName = getName();
    std::string carele = " ";
    std::string test;
    std::string option;
    std::string nompar;
    test = typcham;
    test.resize( 4 );

    if ( test == "ELGA" ) {
        option = "TOU_INI_ELGA";
    } else if ( test == "ELNO" ) {
        option = "TOU_INI_ELNO";
    } else {
        AS_ASSERT( false )
    };
    if ( typcham == test + "_SIEF_R" ) {
        nompar = "PSIEF_R";
    } else if ( typcham == test + "_VARI_R" ) {
        nompar = "PVARI_R";
    } else {
        AS_ASSERT( false )
    };

    if ( carael )
        carele = carael->getName();

    ASTERINTEGER iret = 0;
    auto fed = model->getFiniteElementDescriptor();

    _dofDescription = fed;
    auto dcel = std::make_shared< SimpleFieldOnCells< ASTERINTEGER > >();

    std::string comporName = " ";
    if ( behaviour ) {
        auto compor = behaviour->getBehaviourField();
        comporName = compor->getName();
    }

    CALLO_CESVAR( carele, comporName, fed->getName(), dcel->getName() );
    CALLO_ALCHML( fed->getName(), option, nompar, JeveuxMemoryTypesNames[Permanent], getName(),
                  &iret, dcel->getName() );
    AS_ASSERT( iret == 0 );

    updateValuePointers();
};

template <>
FieldOnCells< ASTERDOUBLE >::FieldOnCells( const ModelPtr &model, const std::string option,
                                           const std::string paraName )
    : FieldOnCells< ASTERDOUBLE >() {
    ASTERINTEGER iret = 0;
    std::string extended = " ";
    auto fed = model->getFiniteElementDescriptor();
    _dofDescription = fed;
    CALLO_ALCHML( fed->getName(), option, paraName, JeveuxMemoryTypesNames[Permanent], getName(),
                  &iret, extended );
    AS_ASSERT( iret == 0 );

    updateValuePointers();
};
