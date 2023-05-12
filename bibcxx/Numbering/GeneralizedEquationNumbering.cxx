/**
 * @file EquationNumbering.cxx
 * @brief Implementation de DOFNumbering
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

#include "Numbering/GeneralizedEquationNumbering.h"

#include "Supervis/ResultNaming.h"

GeneralizedEquationNumbering::GeneralizedEquationNumbering( const std::string name )
    : DataStructure( name, 19, "NUME_EQUA_GENE" ),
      _desc( JeveuxVectorLong( getName() + ".DESC" ) ),
      _nequ( JeveuxVectorLong( getName() + ".NEQU" ) ),
      _refn( JeveuxVectorChar24( getName() + ".REFN" ) ),
      _delg( JeveuxVectorLong( getName() + ".DELG" ) ),
      _orig( JeveuxCollectionLong( getName() + ".ORIG" ) ),
      _componentsOnNodes( getName() + ".PRNO" ),
      _namesOfGroupOfCells( getName() + ".LILI" ),
      _indexationVector( getName() + ".NUEQ" ),
      _nodeAndComponentsIdFromDOF( getName() + ".DEEQ" ) {};

GeneralizedEquationNumbering::GeneralizedEquationNumbering()
    : GeneralizedEquationNumbering( ResultNaming::getNewResultName() ) {};