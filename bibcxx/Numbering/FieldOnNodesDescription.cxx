/**
 * @file FieldOnNodesDescription.cxx
 * @brief Implementation de DOFNumbering
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

#include "Numbering/FieldOnNodesDescription.h"

FieldOnNodesDescription::FieldOnNodesDescription( const std::string name )
    : DataStructure( name, 19, "PROF_CHNO" ), _componentsOnNodes( getName() + ".PRNO" ),
      _namesOfGroupOfCells( getName() + ".LILI" ), _indexationVector( getName() + ".NUEQ" ),
      _nodeAndComponentsNumberFromDOF( getName() + ".DEEQ" ){};

ASTERINTEGER FieldOnNodesDescription::getNumberOfDofs() const {
    return _nodeAndComponentsNumberFromDOF->size() / 2;
};

VectorLong FieldOnNodesDescription::getNodesFromDOF() const {
    const bool retour = _nodeAndComponentsNumberFromDOF->updateValuePointer();
    const ASTERINTEGER nb_eq = this->getNumberOfDofs();

    VectorLong nodes( nb_eq );
    for ( int i_eq = 0; i_eq < nb_eq; i_eq++ )
        nodes[i_eq] = ( *_nodeAndComponentsNumberFromDOF )[2 * i_eq];

    return nodes;
}

/**
 * @brief Mise a jour des pointeurs Jeveux
 * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
 */
bool FieldOnNodesDescription::updateValuePointers() {
    bool retour = _componentsOnNodes->build();
    retour = ( retour && _indexationVector->updateValuePointer() );
    retour = ( retour && _nodeAndComponentsNumberFromDOF->updateValuePointer() );
    return retour;
};

bool FieldOnNodesDescription::operator==( FieldOnNodesDescription &toCompare ) {
    CALL_JEMARQ();
    bool ret = false;

    // TO FIX
    //if ( ( *_componentsOnNodes ) == ( *toCompare._componentsOnNodes ) ) {
        if ( ( *_indexationVector ) == ( *toCompare._indexationVector ) ) {
            if ( ( *_nodeAndComponentsNumberFromDOF ) ==
                 ( *toCompare._nodeAndComponentsNumberFromDOF ) ) {
                if ( ( *_namesOfGroupOfCells ) == ( *toCompare._namesOfGroupOfCells ) ) {
                    ret = true;
                }
            }
        }
    // }
    CALL_JEDEMA();

    return ret;
};
