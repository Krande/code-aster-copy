#ifndef THERMALLOADDESCRIPTION_H_
#define THERMALLOADDESCRIPTION_H_

/**
 * @file ThermalLoad.h
 * @brief Fichier entete de la classe ThermalLoad
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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

#include "DataFields/ConstantFieldOnCells.h"
#include "DataFields/FieldOnCells.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"

/**
 * @class ThermalLoadDescription
 * @brief Classe definissant une charge thermique sd_char_chth
 * @author Jean-Pierre Lefebvre
 */
template < typename ConstantFieldOnCellsType >
class ThermalLoadDescription : public DataStructure {

  public:
    typedef std::shared_ptr< ConstantFieldOnCellsType > ConstantFieldOnCellsTypePtr;

  private:
    /** @brief Modele */
    ModelPtr _model;
    /** @brief Vecteur Jeveux '.MODEL.NOMO' */
    JeveuxVectorChar8 _modelName;
    /** @brief Vecteur Jeveux '.CONVE.VALE' */
    JeveuxVectorChar8 _convection;
    /** @brief Vecteur Jeveux '.EVOL.CHAR' */
    JeveuxVectorChar8 _evolChar;
    /** @brief Vecteur Jeveux '.LIGRE' */
    FiniteElementDescriptorPtr _FEDesc;
    /** @brief Carte '.CIMPO' */
    ConstantFieldOnCellsTypePtr _cimpo;
    /** @brief Carte '.CMULT' */
    ConstantFieldOnCellsRealPtr _cmult;
    /** @brief Carte '.COEFH' */
    ConstantFieldOnCellsTypePtr _coefh;
    /** @brief Carte '.FLUNL' */
    ConstantFieldOnCellsTypePtr _flunl;
    /** @brief Carte '.FLURE' */
    ConstantFieldOnCellsTypePtr _flure;
    /** @brief Carte '.FLUR2' */
    ConstantFieldOnCellsTypePtr _flur2;
    /** @brief Carte '.GRAIN' */
    ConstantFieldOnCellsTypePtr _grain;
    /** @brief Carte '.HECHP' */
    ConstantFieldOnCellsTypePtr _hechp;
    /** @brief Carte '.RAYO' */
    ConstantFieldOnCellsTypePtr _rayo;
    /** @brief Carte '.SOUNL' */
    ConstantFieldOnCellsTypePtr _sounl;
    /** @brief Carte '.SOURE' */
    ConstantFieldOnCellsTypePtr _soure;
    /** @brief Champ '.SOURC' */
    FieldOnCellsRealPtr _sourc;
    /** @brief Carte '.T_EXT' */
    ConstantFieldOnCellsTypePtr _tExt;

  public:
    ThermalLoadDescription( void ) = delete;

    /** @brief Constructeur */
    ThermalLoadDescription( const std::string &name, const ModelPtr &currentModel )
        : DataStructure( name, 13, "CHAR_CHTH" ),
          _model( currentModel ),
          _modelName( getName() + ".MODEL.NOMO" ),
          _convection( getName() + ".CONVE.VALE" ),
          _evolChar( getName() + ".EVOL.CHAR" ),
          _FEDesc( std::make_shared< FiniteElementDescriptor >( getName() + ".LIGRE",
                                                                _model->getMesh() ) ),
          _cimpo( std::make_shared< ConstantFieldOnCellsType >( getName() + ".CIMPO", _FEDesc ) ),
          _cmult( std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CMULT", _FEDesc ) ),
          _coefh( std::make_shared< ConstantFieldOnCellsType >( getName() + ".COEFH", _FEDesc ) ),
          _flunl( std::make_shared< ConstantFieldOnCellsType >( getName() + ".FLUNL", _FEDesc ) ),
          _flure( std::make_shared< ConstantFieldOnCellsType >( getName() + ".FLURE", _FEDesc ) ),
          _flur2( std::make_shared< ConstantFieldOnCellsType >( getName() + ".FLUR2", _FEDesc ) ),
          _grain( std::make_shared< ConstantFieldOnCellsType >( getName() + ".GRAIN", _FEDesc ) ),
          _hechp( std::make_shared< ConstantFieldOnCellsType >( getName() + ".HECHP", _FEDesc ) ),
          _rayo( std::make_shared< ConstantFieldOnCellsType >( getName() + ".RAYO", _FEDesc ) ),
          _sounl( std::make_shared< ConstantFieldOnCellsType >( getName() + ".SOUNL", _FEDesc ) ),
          _soure( std::make_shared< ConstantFieldOnCellsType >( getName() + ".SOURE", _FEDesc ) ),
          _sourc( std::make_shared< FieldOnCellsReal >( getName() + ".SOURC" ) ),
          _tExt( std::make_shared< ConstantFieldOnCellsType >( getName() + ".T_EXT", _FEDesc ) ) {};

    /**
     * @brief Get the finite element descriptor
     */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _FEDesc; };

    ConstantFieldOnCellsTypePtr getImposedField() const { return _cimpo; }

    ConstantFieldOnCellsRealPtr getMultiplicativeField() const { return _cmult; }

    bool hasLoadField( const std::string name ) const {
        if ( name == "COEFH" )
            return ( _coefh && _coefh->exists() );
        else if ( name == "T_EXT" )
            return ( _tExt && _tExt->exists() );
        else if ( name == "FLURE" )
            return ( _flure && _flure->exists() );
        else if ( name == "FLUR2" )
            return ( _flur2 && _flur2->exists() );
        else if ( name == "SOURE" )
            return ( _soure && _soure->exists() );
        else if ( name == "SOURC" )
            return ( _sourc && _sourc->exists() );
        else if ( name == "SOUNL" )
            return _sounl && _sounl->exists();
        else if ( name == "HECHP" )
            return ( _hechp && _hechp->exists() );
        else if ( name == "GRAIN" )
            return ( _grain && _grain->exists() );
        else if ( name == "FLUNL" )
            return _flunl && _flunl->exists();
        else if ( name == "RAYO" )
            return _rayo && _rayo->exists();
        else
            throw std::runtime_error( "Invalid load name : " + name );
    }

    ConstantFieldOnCellsTypePtr getConstantLoadField( const std::string name ) const {
        AS_ASSERT( this->hasLoadField( name ) );
        if ( name == "COEFH" )
            return _coefh;
        else if ( name == "T_EXT" )
            return _tExt;
        else if ( name == "FLURE" )
            return _flure;
        else if ( name == "FLUR2" )
            return _flur2;
        else if ( name == "SOURE" )
            return _soure;
        else if ( name == "SOUNL" )
            return _sounl;
        else if ( name == "HECHP" )
            return _hechp;
        else if ( name == "GRAIN" )
            return _grain;
        else if ( name == "FLUNL" )
            return _flunl;
        else if ( name == "RAYO" )
            return _rayo;
        else
            throw std::runtime_error( "Invalid load name : " + name );
    }

    FieldOnCellsRealPtr getLoadField( const std::string name ) const {
        AS_ASSERT( this->hasLoadField( name ) );
        if ( name == "SOURC" )
            return _sourc;
        else
            throw std::runtime_error( "Invalid load name : " + name );
    }

    bool hasLoadResult() const { return _evolChar->exists(); }

    std::string getLoadResultName() const {
        AS_ASSERT( this->hasLoadResult() );
        _evolChar->updateValuePointer();
        return ( *_evolChar )[0].toString();
    }

    /**
     * @brief Get the model
     */
    ModelPtr getModel() const { return _model; };

    /**
     * @brief Get the mesh
     */
    BaseMeshPtr getMesh() const { return _model->getMesh(); };

    bool build() {
        _FEDesc->build();

        return true;
    };
};

/**********************************************************
 *  Explicit instantiation of template classes
 **********************************************************/

/** @typedef ThermalLoadDescriptionReal Class d'une charge mécanique réelle */
typedef ThermalLoadDescription< ConstantFieldOnCellsReal > ThermalLoadDescriptionReal;
/** @typedef ThermalLoadDescriptionFunc Class d'une charge mécanique de fonctions */
typedef ThermalLoadDescription< ConstantFieldOnCellsChar24 > ThermalLoadDescriptionFunction;

template < typename ConstantFieldOnCellsType >
using ThermalLoadDescriptionPtr =
    std::shared_ptr< ThermalLoadDescription< ConstantFieldOnCellsType > >;

#endif /* THERMALLOADDESCRIPTION_H_ */
