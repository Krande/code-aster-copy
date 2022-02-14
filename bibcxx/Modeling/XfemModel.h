#ifndef XFEMMODEL_H_
#define XFEMMODEL_H_

/**
 * @file XfemModel.h
 * @brief Header for class XfemModel
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

#include "astercxx.h"

#include "DataFields/DataField.h"

/** @class XfemField */
class XfemField : public DataField {
  private:
    /** @brief Vecteur Jeveux '.CELD' */
    JeveuxVectorLong _celd;
    /** @brief Vecteur Jeveux '.CELK' */
    JeveuxVectorChar24 _celk;
    /** @brief Vecteur Jeveux '.CELV' */
    JeveuxVector< ASTERDOUBLE > _celv;

  public:
    /** @typedef XfemFieldOnCellsPtr */
    typedef std::shared_ptr< XfemField > XfemFieldPtr;

    /** @brief Constructor */
    XfemField( const std::string name )
        : DataField( name, "CHAM_ELEM" ),
          _celd( JeveuxVectorLong( getName() + ".CELD" ) ),
          _celk( JeveuxVectorChar24( getName() + ".CELK" ) ),
          _celv( JeveuxVector< ASTERDOUBLE >( getName() + ".CELV" ) ){};

    /** @brief Constructor */
    XfemField() : XfemField( ResultNaming::getNewResultName() ){};

    /** @brief Deallocate JEVEUX pointers */
    void deallocate() {
        _celd->deallocate();
        _celk->deallocate();
        _celv->deallocate();
    };
};

typedef std::shared_ptr< XfemField > XfemFieldPtr;

/**
 * @class XfemModel
 * @brief Datastructure for XfemModel (AFFE_MODELE)
 */
class XfemModel {
  private:
    typedef std::shared_ptr< XfemField > XfemFieldPtr;
    typedef std::map< std::string, XfemFieldPtr > listFields;
    /** Fields for XFEM */
    listFields _listfields;

  public:
    XfemModel( const std::string modelName ) { createFields( modelName ); };

    void createFields( const std::string modelName ) {
        XfemFieldPtr pinttoPtr = std::make_shared< XfemField >( modelName + ".TOPOSE.PIN" );
        _listfields.insert( listFields::value_type( "PINTTO", pinttoPtr ) );
        XfemFieldPtr cnsetoPtr = std::make_shared< XfemField >( modelName + ".TOPOSE.CNS" );
        _listfields.insert( listFields::value_type( "CNSETO", cnsetoPtr ) );
        XfemFieldPtr heavtoPtr = std::make_shared< XfemField >( modelName + ".TOPOSE.HEA" );
        _listfields.insert( listFields::value_type( "HEAVTO", heavtoPtr ) );
        XfemFieldPtr lonchaPtr = std::make_shared< XfemField >( modelName + ".TOPOSE.LON" );
        _listfields.insert( listFields::value_type( "LONCHA", lonchaPtr ) );
        XfemFieldPtr baslocPtr = std::make_shared< XfemField >( modelName + ".BASLOC" );
        _listfields.insert( listFields::value_type( "BASLOC", baslocPtr ) );
        XfemFieldPtr lsnPtr = std::make_shared< XfemField >( modelName + ".LNNO" );
        _listfields.insert( listFields::value_type( "LSN", lsnPtr ) );
        XfemFieldPtr lstPtr = std::make_shared< XfemField >( modelName + ".LTNO" );
        _listfields.insert( listFields::value_type( "LST", lstPtr ) );
        XfemFieldPtr stanoPtr = std::make_shared< XfemField >( modelName + ".STNO" );
        _listfields.insert( listFields::value_type( "STANO", stanoPtr ) );
        XfemFieldPtr pmiltoPtr = std::make_shared< XfemField >( modelName + ".TOPOSE.PMI" );
        _listfields.insert( listFields::value_type( "PMILT", pmiltoPtr ) );
        XfemFieldPtr fissnoPtr = std::make_shared< XfemField >( modelName + ".FISSNO" );
        _listfields.insert( listFields::value_type( "FISSNO", fissnoPtr ) );
        XfemFieldPtr heavnoPtr = std::make_shared< XfemField >( modelName + ".TOPONO.HNO" );
        _listfields.insert( listFields::value_type( "HEAVNO", heavnoPtr ) );
        XfemFieldPtr heavfaPtr = std::make_shared< XfemField >( modelName + ".TOPONO.HFA" );
        _listfields.insert( listFields::value_type( "HEAVFA", heavfaPtr ) );
    };

    XfemFieldPtr getField( const std::string fieldType ) const {
        return _listfields.at( fieldType );
    };
};

typedef std::shared_ptr< XfemModel > XfemModelPtr;

#endif /* XFEMMODEL_H_ */
