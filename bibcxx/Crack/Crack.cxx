/**
 * @file Crack.cxx
 * @brief Implementation de Crack
 * @author Nicolas Pignet
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
      _ltno( std::make_shared< FieldOnNodesReal >( getName() + ".LTNO      " ) ),
      _lnno( std::make_shared< FieldOnNodesReal >( getName() + ".LNNO      " ) ),
      _basLoc( std::make_shared< FieldOnNodesReal >( getName() + ".BASLOC    " ) ),
      _basNof( JeveuxVectorReal( getName() + ".BASNOF" ) ),
      _coorfond( JeveuxVectorReal( getName() + ".COORFOND" ) ),
      _absfon( JeveuxVectorReal( getName() + ".ABSFON" ) ) {
    _ltno->setDescription( std::make_shared< EquationNumbering >( getName() + ".LTNO.NUMEQ" ) );
    _lnno->setDescription( _ltno->getDescription() );
    _basLoc->setDescription( std::make_shared< EquationNumbering >( getName() + ".BASL.NUMEQ" ) );
};

void Crack::updateValuePointers() {
    _info->updateValuePointer();
    _basNof->updateValuePointer();
    _coorfond->updateValuePointer();
}

std::string Crack::getCrackTipCellsType() {
    this->updateValuePointers();
    return strip( ( *_info )[4].toString() );
}

std::string Crack::getUpperLipGroupName() {
    this->updateValuePointers();
    return strip( ( *_info )[5].toString() );
}

std::string Crack::getLowerLipGroupName() {
    this->updateValuePointers();
    return strip( ( *_info )[6].toString() );
}

const JeveuxVectorReal Crack::getCrackFrontBasis() {
    this->updateValuePointers();
    return _basNof;
};

const VectorString Crack::getCrackFrontNodes() {
    VectorString crackFrontNodes;
    _fondNoeu->updateValuePointer();
    crackFrontNodes.reserve( _fondNoeu->size() );
    for ( auto &node : _fondNoeu )
        crackFrontNodes.push_back( strip( node ) );
    return crackFrontNodes;
};

const VectorString Crack::getLowerNormNodes() {
    VectorString lowerNormNodes;
    _infNormNoeud->updateValuePointer();
    lowerNormNodes.reserve( _infNormNoeud->size() );
    for ( auto &node : _infNormNoeud )
        lowerNormNodes.push_back( strip( node ) );
    return lowerNormNodes;
};

const VectorString Crack::getUpperNormNodes() {
    VectorString upperNormNodes;
    _supNormNoeu->updateValuePointer();
    upperNormNodes.reserve( _supNormNoeu->size() );
    for ( auto &node : _supNormNoeu )
        upperNormNodes.push_back( strip( node ) );
    return upperNormNodes;
};

const VectorString Crack::getLowerNormNodes2() {
    VectorString lowerNormNodes;
    _infNormNoeud2->updateValuePointer();
    lowerNormNodes.reserve( _infNormNoeud2->size() );
    for ( auto &node : _infNormNoeud2 )
        lowerNormNodes.push_back( strip( node ) );
    return lowerNormNodes;
};

const VectorString Crack::getUpperNormNodes2() {
    VectorString upperNormNodes;
    _supNormNoeu2->updateValuePointer();
    upperNormNodes.reserve( _supNormNoeu2->size() );
    for ( auto &node : _supNormNoeu2 )
        upperNormNodes.push_back( strip( node ) );
    return upperNormNodes;
};

const JeveuxVectorReal Crack::getCrackFrontPosition() {
    this->updateValuePointers();
    return _coorfond;
};

const JeveuxVectorReal Crack::getCrackFrontAbsCurv() {
    this->updateValuePointers();
    return _absfon;
};

const FieldOnNodesRealPtr Crack::getCrackFrontNodeBasis() { return _basLoc; };

const JeveuxVectorReal Crack::getCrackFrontRadius() { return _fondTailleR; };

const JeveuxVectorReal Crack::getNormal() { return _normale; };

bool Crack::isSymmetric() {
    this->updateValuePointers();
    return strip( ( *_info )[0].toString() ) == "OUI";
};

std::string Crack::getConfigInit() {
    this->updateValuePointers();
    return strip( ( *_info )[1].toString() );
};
