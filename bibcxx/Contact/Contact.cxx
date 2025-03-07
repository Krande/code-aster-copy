/**
 * @file Contact.cxx
 * @brief Implementation de Contact
 * @author Nicolas Sellenet
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

#include "Contact/Contact.h"

Contact::Contact( const std::string name, const ModelPtr model )
    : DataStructure( name, 8, "CHAR_CONTACT" ),
      _model( model ),
      _FEDesc( std::make_shared< FiniteElementDescriptor >( getName() + ".CONT.LIGRE",
                                                            _model->getMesh() ) ),
      _isEmpty( false ),
      _model_name( JeveuxVectorChar8( getName() + ".CHME.MODEL.NOMO" ) ),
      _integer_params( JeveuxVectorLong( getName() + ".PARACI" ) ),
      _real_params( JeveuxVectorReal( getName() + ".PARACR" ) ),
      _contact_type( JeveuxVectorChar8( getName() + ".TYPE" ) ),
      _ndim_unilate( JeveuxVectorLong( getName() + ".UNILATE.NDIMCU" ) ),
      _cmpg_unilate( JeveuxVectorChar8( getName() + ".UNILATE.CMPGCU" ) ),
      _coed_unilate( JeveuxVectorChar8( getName() + ".UNILATE.COED" ) ),
      _coeg_unilate( JeveuxVectorChar8( getName() + ".UNILATE.COEG" ) ),
      _lisnoe_unilate( JeveuxVectorLong( getName() + ".UNILATE.LISNOE" ) ),
      _poinoe_unilate( JeveuxVectorLong( getName() + ".UNILATE.POINOE" ) ),
      _coefpe_unilate( JeveuxVectorReal( getName() + ".UNILATE.COEFPE" ) ),
      _ndimco( JeveuxVectorLong( getName() + ".CONTACT.NDIMCO" ) ),
      _methco( JeveuxVectorLong( getName() + ".CONTACT.METHCO" ) ),
      _dirapp( JeveuxVectorChar8( getName() + ".CONTACT.DIRAPP" ) ),
      _dirnor( JeveuxVectorChar8( getName() + ".CONTACT.DIRNOR" ) ),
      _jfo1co( JeveuxVectorChar8( getName() + ".CONTACT.JFO1CO" ) ),
      _jfo2co( JeveuxVectorChar8( getName() + ".CONTACT.JFO2CO" ) ),
      _toleco( JeveuxVectorReal( getName() + ".CONTACT.TOLECO" ) ),
      _jeucoq( JeveuxVectorReal( getName() + ".CONTACT.JEUCOQ" ) ),
      _jeupou( JeveuxVectorReal( getName() + ".CONTACT.JEUPOU" ) ),
      _pzoneco( JeveuxVectorLong( getName() + ".CONTACT.PZONECO" ) ),
      _psumaco( JeveuxVectorLong( getName() + ".CONTACT.PSUMACO" ) ),
      _psunoco( JeveuxVectorLong( getName() + ".CONTACT.PSUNOCO" ) ),
      _mailco( JeveuxVectorLong( getName() + ".CONTACT.MAILCO" ) ),
      _noeuco( JeveuxVectorLong( getName() + ".CONTACT.NOEUCO" ) ),
      _manoco( JeveuxVectorLong( getName() + ".CONTACT.MANOCO" ) ),
      _pmanoco( JeveuxVectorLong( getName() + ".CONTACT.PMANOCO" ) ),
      _nomaco( JeveuxVectorLong( getName() + ".CONTACT.NOMACO" ) ),
      _pnomaco( JeveuxVectorLong( getName() + ".CONTACT.PNOMACO" ) ),
      _pssnoco( JeveuxVectorLong( getName() + ".CONTACT.PSSNOCO" ) ),
      _ssnoco( JeveuxVectorLong( getName() + ".CONTACT.SSNOCO" ) ),
      _typeno( JeveuxVectorLong( getName() + ".CONTACT.TYPENO" ) ),
      _typema( JeveuxVectorLong( getName() + ".CONTACT.TYPEMA" ) ),
      _maescl( JeveuxVectorLong( getName() + ".CONTACT.MAESCL" ) ),
      _caradf( JeveuxVectorReal( getName() + ".CONTACT.CARADF" ) ),
      _caracf( JeveuxVectorReal( getName() + ".CONTACT.CARACF" ) ),
      _psanofr( JeveuxVectorLong( getName() + ".CONTACT.PSANOFR" ) ),
      _sanofr( JeveuxVectorLong( getName() + ".CONTACT.SANOFR" ) ),
      _exclfr( JeveuxVectorLong( getName() + ".CONTACT.EXCLFR" ) ),
      _caraxf( JeveuxVectorReal( getName() + ".CONTACT.CARAXF" ) ),
      _xfimai( JeveuxVectorChar8( getName() + ".CONTACT.XFIMAI" ) ),
      _xnrell( JeveuxVectorChar24( getName() + ".CONTACT.XNRELL" ) ),
      _maescx( JeveuxVectorLong( getName() + ".CONTACT.MAESCX" ) ),
      _ptrdclc( JeveuxVectorLong( getName() + ".CONTACT.PTRDCLC" ) ) {};

Contact::Contact( const ModelPtr model ) : Contact( ResultNaming::getNewResultName(), model ) {};

ModelPtr Contact::getModel() const { return _model; };

FiniteElementDescriptorPtr Contact::getFiniteElementDescriptor() const { return _FEDesc; };
