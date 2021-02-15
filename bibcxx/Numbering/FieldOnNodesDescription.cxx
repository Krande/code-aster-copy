/**
 * @file FieldOnNodesDescription.cxx
 * @brief Implementation de DOFNumbering
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

FieldOnNodesDescriptionClass::FieldOnNodesDescriptionClass( const JeveuxMemory memType )
    : DataStructure( ResultNaming::getNewResultName(), 19, "PROF_CHNO", memType ),
      _componentsOnNodes( getName() + ".PRNO" ), _namesOfGroupOfCells( getName() + ".LILI" ),
      _indexationVector( getName() + ".NUEQ" ),
      _nodeAndComponentsNumberFromDOF( getName() + ".DEEQ" ){};

FieldOnNodesDescriptionClass::FieldOnNodesDescriptionClass( const std::string name,
                                                            const JeveuxMemory memType )
    : DataStructure( name, 19, "PROF_CHNO", memType ), _componentsOnNodes( getName() + ".PRNO" ),
      _namesOfGroupOfCells( getName() + ".LILI" ), _indexationVector( getName() + ".NUEQ" ),
      _nodeAndComponentsNumberFromDOF( getName() + ".DEEQ" ){};



ASTERINTEGER FieldOnNodesDescriptionClass::getNumberOfDofs() const
{
    return _nodeAndComponentsNumberFromDOF->size() / 2;
};

VectorLong FieldOnNodesDescriptionClass::getNodesFromDOF() const
{
    const bool retour = _nodeAndComponentsNumberFromDOF->updateValuePointer();
    const ASTERINTEGER nb_eq = this->getNumberOfDofs();

    VectorLong nodes(nb_eq);
    for(int i_eq = 0; i_eq < nb_eq; i_eq++)
        nodes[i_eq] = (*_nodeAndComponentsNumberFromDOF)[2 * i_eq];

    return nodes;
}
