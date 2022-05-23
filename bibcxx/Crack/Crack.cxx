/**
 * @file Crack.cxx
 * @brief Implementation de Crack
 * @author Nicolas Pignet
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

/* person_in_charge: nicolas.pignet at edf.fr */

#include "astercxx.h"

#include "Crack/Crack.h"

Crack::Crack( const std::string name )
    : DataStructure( name, 8, "FOND_FISSURE" ),
      _levreInfMail( JeveuxVectorChar8( getName() + ".LEVREINF.MAIL" ) ),
      _normale( JeveuxVectorReal( getName() + ".NORMALE" ) ),
      _fondNoeu( JeveuxVectorChar8( getName() + ".FOND.NOEU" ) ),
      _infNormNoeud( JeveuxVectorChar8( getName() + ".INFNORM.NOEU" ) ),
      _supNormNoeu( JeveuxVectorChar8( getName() + ".SUPNORM.NOEU" ) ),
      _infNormNoeud2( JeveuxVectorChar8( getName() + ".INFNORM.NOEU2" ) ),
      _supNormNoeu2( JeveuxVectorChar8( getName() + ".SUPNORM.NOEU2" ) ),
      _levreSupMail( JeveuxVectorChar8( getName() + ".LEVRESUP.MAIL" ) ),
      _info( JeveuxVectorChar24( getName() + ".INFO" ) ),
      _fondTailleR( JeveuxVectorReal( getName() + ".FOND.TAILLE_R" ) ),
      _abscur( JeveuxVectorReal( getName() + ".ABSCUR" ) ),
      _ltno( new FieldOnNodesReal( getName() + ".LTNO      " ) ),
      _lnno( new FieldOnNodesReal( getName() + ".LNNO      " ) ),
      _basLoc( new FieldOnNodesReal( getName() + ".BASLOC    " ) ),
      _basNof( JeveuxVectorReal( getName() + ".BASNOF" ) ),
      _absfon( JeveuxVectorReal( getName() + ".ABSFON" ) ){};

void Crack::updateValuePointers( ){
    _info->updateValuePointer();
}

std::string Crack::getCrackTipCellsType( ){
  this->updateValuePointers();
  return trim( ( *_info )[4].toString() );
}

std::string Crack::getUpperLipGroupName( ){
  this->updateValuePointers();
  return trim( ( *_info )[5].toString() );
}

std::string Crack::getLowerLipGroupName( ){
  this->updateValuePointers();
  return trim( ( *_info )[6].toString() );
}
