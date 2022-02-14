/**
 * @file ExternalStateVariables.cxx
 * @brief Implementation de ExternalStateVariable
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

#include "Materials/ExternalStateVariablesBuilder.h"

#include "Materials/BaseExternalStateVariables.h"

ExternalStateVariablesBuilder::ExternalStateVariablesBuilder()
    : DataStructure( "VARI_COM", 14 ),
      _varcVale( new FieldOnCellsReal( getName() + ".TOUT" ) ),
      _timeField( new ConstantFieldOnCellsReal( getName() + ".INST", _model->getMesh() ) ),
      _time( -1.0 ),
      _pTot( false ),
      _hydr( false ),
      _sech( false ),
      _temp( false ){};

ExternalStateVariablesBuilder::ExternalStateVariablesBuilder(
    const ModelPtr &currModel, const MaterialFieldPtr &currMaterialField,
    const ElementaryCharacteristicsPtr &currElemChara, const CodedMaterialPtr &currCodedMaterial )
    : DataStructure( "VARI_COM", 14 ),
      _model( currModel ),
      _materialField( currMaterialField ),
      _codedMaterial( currCodedMaterial ),
      _elemChara( currElemChara ),
      _varcRefe( new FieldOnCellsReal( _model->getName() + ".CHVCREF" ) ),
      _varcVale( new FieldOnCellsReal( getName() + ".TOUT" ) ),
      _timeField( new ConstantFieldOnCellsReal( getName() + ".INST", _model->getMesh() ) ),
      _time( -1.0 ),
      _pTot( _materialField->hasExternalStateVariables( "PTOT" ) ),
      _hydr( _materialField->hasExternalStateVariables( "HYDR" ) ),
      _sech( _materialField->hasExternalStateVariables( "SECH" ) ),
      _temp( _materialField->hasExternalStateVariables( "TEMP" ) ) {
    std::string modelName( _model->getName(), 0, 8 ),
        materialFieldName( _materialField->getName(), 0, 8 );
    std::string elemCharaName( ' ', 8 );
    if ( _elemChara != nullptr )
        elemCharaName = std::string( _elemChara->getName(), 0, 8 );
    CALLO_VRCREF( modelName, materialFieldName, elemCharaName, _varcRefe->getName() );
};

void ExternalStateVariablesBuilder::build( const ASTERDOUBLE &time ) {
    _time = time;
    _varcVale->deallocate();
    _timeField->deallocate();

    std::string modelName( _model->getName(), 0, 8 );
    std::string materialFieldName( _materialField->getName(), 0, 8 );
    std::string elemCharaName( 8, ' ' );
    if ( _elemChara != nullptr )
        elemCharaName = std::string( _elemChara->getName(), 0, 8 );
    std::string out( ' ', 2 );
    CALLO_VRCINS_WRAP( modelName, materialFieldName, elemCharaName, &time, _varcVale->getName(),
                       out );

    std::string comp( "INST_R" );
    _timeField->allocate( comp );
    ConstantFieldOnZone a( _model->getMesh() );
    ConstantFieldValues< ASTERDOUBLE > b( { "INST" }, { time } );
    _timeField->setValueOnZone( a, b );
};

FieldOnNodesRealPtr ExternalStateVariablesBuilder::computeExternalStateVariablesLoad(
    const BaseDOFNumberingPtr &dofNUM ) {
    const auto &codedMaterial = _codedMaterial->getCodedMaterialField();
    std::string modelName( _model->getName(), 0, 8 );
    std::string elemCharaName( 8, ' ' );
    if ( _elemChara != nullptr )
        elemCharaName = std::string( _elemChara->getName(), 0, 8 );

    const auto compor = _materialField->getBehaviourField();
    std::string numName( dofNUM->getName() ), out( 24, ' ' );
    ASTERINTEGER a = 0, b = 0, c = 0, d = 0;
    if ( _hydr )
        a = 1;
    if ( _sech )
        b = 1;
    if ( _temp )
        c = 1;
    if ( _pTot )
        d = 1;
    CALLO_CACHVC( modelName, _materialField->getName(), codedMaterial->getName(), elemCharaName,
                  numName, compor->getName(), getName(), out, &a, &b, &c, &d );
    JeveuxVectorChar24 vectOut( out );
    vectOut->updateValuePointer();
    FieldOnNodesRealPtr toReturn( new FieldOnNodesReal( ( *vectOut )[0].toString() ) );
    toReturn->setDescription( dofNUM->getDescription() );
    toReturn->setMesh( dofNUM->getMesh() );
    toReturn->build();
    toReturn->updateValuePointers();

    return toReturn;
};
