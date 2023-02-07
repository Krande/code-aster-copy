#ifndef MECHANICALLOADDESCRIPTION_H_
#define MECHANICALLOADDESCRIPTION_H_

/**
 * @file MechanicalLoadDescription.h
 * @author Natacha Bereux
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

#include "astercxx.h"

#include "DataFields/ConstantFieldOnCells.h"
#include "DataStructures/DataStructure.h"
#include "Meshes/BaseMesh.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"

#include <string>

/**
 * @class MechanicalLoadDescription
 * @brief Encapsulation of sd_char_chme
 */
template < typename ConstantFieldOnCellsType >
class MechanicalLoadDescription : public DataStructure {

  public:
    typedef std::shared_ptr< ConstantFieldOnCellsType > ConstantFieldOnCellsTypePtr;

  private:
    /** @brief Modele */
    ModelPtr _model;
    /** @brief Vecteur Jeveux '.TEMPE.TEMP' */
    JeveuxVectorChar8 _temperatureField;
    /** @brief Vecteur Jeveux '.MODEL.NOMO' */
    JeveuxVectorChar8 _modelName;
    /** @brief Vecteur Jeveux '.VEASS' */
    JeveuxVectorChar8 _nameOfAssemblyVector;
    /** @brief Vecteur Jeveux '.VEISS' */
    JeveuxVectorChar8 _veiss;
    /** @brief Vecteur Jeveux '.EVOL.CHAR' */
    JeveuxVectorChar8 _evolChar;
    /** @brief Vecteur Jeveux '.LIGRE' */
    FiniteElementDescriptorPtr _FEDesc;
    /** @brief Carte '.CIMPO' */
    ConstantFieldOnCellsTypePtr _cimpo;
    /** @brief Carte '.CMULT' */
    ConstantFieldOnCellsRealPtr _cmult;
    /** @brief Carte '.DPGEN' */
    ConstantFieldOnCellsTypePtr _dpgen;
    /** @brief Carte '.EPSIN' */
    ConstantFieldOnCellsTypePtr _epsin;
    /** @brief Carte '.F1D1D' */
    ConstantFieldOnCellsTypePtr _f1d1d;
    /** @brief Carte '.F1D2D' */
    ConstantFieldOnCellsTypePtr _f1d2d;
    /** @brief Carte '.F1D3D' */
    ConstantFieldOnCellsTypePtr _f1d3d;
    /** @brief Carte '.F2D2D' */
    ConstantFieldOnCellsTypePtr _f2d2d;
    /** @brief Carte '.F2D3D' */
    ConstantFieldOnCellsTypePtr _f2d3d;
    /** @brief Carte '.F3D3D' */
    ConstantFieldOnCellsTypePtr _f3d3d;
    /** @brief Carte '.FCO2D' */
    ConstantFieldOnCellsTypePtr _fco2d;
    /** @brief Carte '.FCO3D' */
    ConstantFieldOnCellsTypePtr _fco3d;
    /** @brief Carte '.FELEC' */
    ConstantFieldOnCellsTypePtr _felec;
    /** @brief Carte '.FL101' */
    ConstantFieldOnCellsTypePtr _fl101;
    /** @brief Carte '.FL102' */
    ConstantFieldOnCellsTypePtr _fl102;
    /** @brief Carte '.FORNO' */
    ConstantFieldOnCellsTypePtr _forno;
    /** @brief Carte '.IMPED' */
    ConstantFieldOnCellsTypePtr _imped;
    /** @brief Carte '.PESAN' */
    ConstantFieldOnCellsTypePtr _pesan;
    /** @brief Carte '.PRESS' */
    ConstantFieldOnCellsTypePtr _press;
    /** @brief Carte '.ROTAT' */
    ConstantFieldOnCellsTypePtr _rotat;
    /** @brief Carte '.SIGIN' */
    ConstantFieldOnCellsTypePtr _sigin;
    /** @brief Carte '.SIINT' */
    ConstantFieldOnCellsTypePtr _siint;
    /** @brief Carte '.VNOR' */
    ConstantFieldOnCellsTypePtr _vnor;
    /** @brief Carte '.ONDPL' */
    ConstantFieldOnCellsTypePtr _ondpl;
    /** @brief Carte '.ONDPR' */
    ConstantFieldOnCellsTypePtr _ondpr;
    /** @brief Carte '.ONDE' */
    ConstantFieldOnCellsTypePtr _onde;

  public:
    MechanicalLoadDescription( void ) = delete;

    /**
     * @brief Constructor
     */
    MechanicalLoadDescription( const std::string name, const ModelPtr &currentModel )
        : DataStructure( name, 13, "CHAR_CHME" ),
          _model( currentModel ),
          _temperatureField( getName() + ".TEMPE.TEMP" ),
          _modelName( getName() + ".MODEL.NOMO" ),
          _nameOfAssemblyVector( getName() + ".VEASS" ),
          _veiss( getName() + ".VEISS" ),
          _evolChar( getName() + ".EVOL.CHAR" ),
          _FEDesc( std::make_shared< FiniteElementDescriptor >( getName() + ".LIGRE",
                                                                _model->getMesh() ) ),
          _cimpo( std::make_shared< ConstantFieldOnCellsType >( getName() + ".CIMPO", _FEDesc ) ),
          _cmult( std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CMULT", _FEDesc ) ),
          _dpgen( std::make_shared< ConstantFieldOnCellsType >( getName() + ".DPGEN", _FEDesc ) ),
          _epsin( std::make_shared< ConstantFieldOnCellsType >( getName() + ".EPSIN", _FEDesc ) ),
          _f1d1d( std::make_shared< ConstantFieldOnCellsType >( getName() + ".F1D1D", _FEDesc ) ),
          _f1d2d( std::make_shared< ConstantFieldOnCellsType >( getName() + ".F1D2D", _FEDesc ) ),
          _f2d2d( std::make_shared< ConstantFieldOnCellsType >( getName() + ".F2D2D", _FEDesc ) ),
          _f1d3d( std::make_shared< ConstantFieldOnCellsType >( getName() + ".F1D3D", _FEDesc ) ),
          _f2d3d( std::make_shared< ConstantFieldOnCellsType >( getName() + ".F2D3D", _FEDesc ) ),
          _f3d3d( std::make_shared< ConstantFieldOnCellsType >( getName() + ".F3D3D", _FEDesc ) ),
          _fco2d( std::make_shared< ConstantFieldOnCellsType >( getName() + ".FCO2D", _FEDesc ) ),
          _fco3d( std::make_shared< ConstantFieldOnCellsType >( getName() + ".FCO3D", _FEDesc ) ),
          _felec( std::make_shared< ConstantFieldOnCellsType >( getName() + ".FELEC", _FEDesc ) ),
          _fl101( std::make_shared< ConstantFieldOnCellsType >( getName() + ".FL101", _FEDesc ) ),
          _fl102( std::make_shared< ConstantFieldOnCellsType >( getName() + ".FL102", _FEDesc ) ),
          _forno( std::make_shared< ConstantFieldOnCellsType >( getName() + ".FORNO", _FEDesc ) ),
          _imped( std::make_shared< ConstantFieldOnCellsType >( getName() + ".IMPED", _FEDesc ) ),
          _pesan( std::make_shared< ConstantFieldOnCellsType >( getName() + ".PESAN", _FEDesc ) ),
          _press( std::make_shared< ConstantFieldOnCellsType >( getName() + ".PRESS", _FEDesc ) ),
          _rotat( std::make_shared< ConstantFieldOnCellsType >( getName() + ".ROTAT", _FEDesc ) ),
          _sigin( std::make_shared< ConstantFieldOnCellsType >( getName() + ".SIGIN", _FEDesc ) ),
          _siint( std::make_shared< ConstantFieldOnCellsType >( getName() + ".SIINT", _FEDesc ) ),
          _vnor( std::make_shared< ConstantFieldOnCellsType >( getName() + ".VNOR", _FEDesc ) ),
          _onde( std::make_shared< ConstantFieldOnCellsType >( getName() + ".ONDE", _FEDesc ) ),
          _ondpl( std::make_shared< ConstantFieldOnCellsType >( getName() + ".ONDPL", _FEDesc ) ),
          _ondpr( std::make_shared< ConstantFieldOnCellsType >( getName() + ".ONDPR", _FEDesc ) ){};

    /**
     * @brief Get the finite element descriptor
     */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _FEDesc; };

    BaseMeshPtr getMesh() const { return _model->getMesh(); }

    ModelPtr getModel() const { return _model; }

    ConstantFieldOnCellsTypePtr getImposedField() const { return _cimpo; }
    ConstantFieldOnCellsRealPtr getMultiplicativeField() const { return _cmult; }

    bool hasLoadField( const std::string load_name ) const {
        if ( load_name == "IMPED" ) {
            return ( _imped && _imped->exists() );
        } else if ( load_name == "ROTAT" ) {
            return ( _rotat && _rotat->exists() );
        } else if ( load_name == "ONDE" ) {
            return ( _onde && _onde->exists() );
        } else {
            AS_ASSERT( false );
        }

        return false;
    }

    ConstantFieldOnCellsTypePtr getConstantLoadField( const std::string name ) const {
        if ( name == "IMPED" ) {
            return _imped;
        } else if ( name == "ROTAT" ) {
            return _rotat;
        } else if ( name == "ONDE" ) {
            return _onde;
        } else {
            AS_ASSERT( false );
        }

        return nullptr;
    }

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return true si la mise a jour s'est bien deroulee, false sinon
     */
    void updateValuePointers() {
        _temperatureField->updateValuePointer();
        _modelName->updateValuePointer();
        _nameOfAssemblyVector->updateValuePointer();
        _veiss->updateValuePointer();
        _evolChar->updateValuePointer();
        _cimpo->updateValuePointers();
        _cmult->updateValuePointers();
        _dpgen->updateValuePointers();
        _epsin->updateValuePointers();
        _f1d2d->updateValuePointers();
        _f1d3d->updateValuePointers();
        _f2d3d->updateValuePointers();
        _fco2d->updateValuePointers();
        _fco3d->updateValuePointers();
        _felec->updateValuePointers();
        _fl101->updateValuePointers();
        _fl102->updateValuePointers();
        _forno->updateValuePointers();
        _imped->updateValuePointers();
        _pesan->updateValuePointers();
        _press->updateValuePointers();
        _sigin->updateValuePointers();
        _rotat->updateValuePointers();
        _siint->updateValuePointers();
        _vnor->updateValuePointers();
        _ondpl->updateValuePointers();
        _ondpr->updateValuePointers();
    };

    bool build() {
        _FEDesc->build();

        return true;
    };
};

/**********************************************************/
/*  Explicit instantiation of template classes
/**********************************************************/

/** @typedef MechanicalLoadDescriptionReal Class d'une charge mécanique réelle */
typedef MechanicalLoadDescription< ConstantFieldOnCellsReal > MechanicalLoadDescriptionReal;
/** @typedef MechanicalLoadDescriptionFunc Class d'une charge mécanique de fonctions */
typedef MechanicalLoadDescription< ConstantFieldOnCellsChar24 > MechanicalLoadDescriptionFunction;
/** @typedef MechanicalLoadDescriptionComplex Class d'une charge mécanique de complexe */
typedef MechanicalLoadDescription< ConstantFieldOnCellsComplex > MechanicalLoadDescriptionComplex;

template < typename ConstantFieldOnCellsType >
using MechanicalLoadDescriptionPtr =
    std::shared_ptr< MechanicalLoadDescription< ConstantFieldOnCellsType > >;

#endif /* MECHANICALLOADDESCRIPTION_H_ */
