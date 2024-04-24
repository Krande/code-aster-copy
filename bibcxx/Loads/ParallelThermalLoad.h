
#include "astercxx.h"

#ifdef ASTER_HAVE_MPI

#ifndef PARALLELTHERMALLOAD_H_
#define PARALLELTHERMALLOAD_H_

/**
 * @file ParallelThermalLoad.h
 * @brief Fichier entete de la classe ParallelThermalLoad
 * @author Nicolas Sellenet
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

#include "DataFields/ConstantFieldOnCells.h"
#include "DataStructures/DataStructure.h"
#include "Loads/ThermalLoad.h"
#include "Meshes/MeshExplorer.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Modeling/ParallelFiniteElementDescriptor.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Supervis/ResultNaming.h"
#include "Utilities/SyntaxSaver.h"

/**
 * @class ParallelThermalLoad
 * @brief Classe definissant une charge dualisée parallèle
 */
template < typename ConstantFieldOnCellsType >
class ParallelThermalLoad : public DataStructure {
  public:
    using ConstantFieldOnCellsTypePtr = std::shared_ptr< ConstantFieldOnCellsType >;

    using ParallelThermalLoadPtr = std::shared_ptr< ParallelThermalLoad >;

  private:
    template < typename ConstantFieldOnCellsType2Ptr >
    void transferConstantFieldOnCells( const ConstantFieldOnCellsType2Ptr &fieldIn,
                                       ConstantFieldOnCellsType2Ptr &fieldOut ) {
        const auto &toKeep = _FEDesc->getVirtualCellsToKeep();

        std::string savedName( "" );
        fieldIn->build();
        fieldOut->allocate( fieldIn );
        const auto sizeFieldIn = ( *fieldIn ).size();
        const auto vect_resu = fieldIn->getValues();

        fieldIn->updateValuePointers();
        for ( int pos = 0; pos < sizeFieldIn; ++pos ) {
            const auto &zone = fieldIn->getZoneDescription( pos );
            const auto &curFEDesc = zone.getFiniteElementDescriptor();
            if ( curFEDesc->getName() != savedName && savedName != "" ) {
                std::string a( "Different FiniteElementDescriptor in one ConstantFieldOnCells is "
                               "not allowed" );
                throw std::runtime_error( a );
            }
            savedName = curFEDesc->getName();

            const auto &listCells = zone.getListOfCells();
            VectorLong toCopy;
            toCopy.reserve( listCells.size() );
            for ( const auto &num : listCells ) {
                if ( toKeep[-num - 1] != 1 )
                    toCopy.push_back( toKeep[-num - 1] );
            }

            if ( toCopy.size() != 0 ) {
                const auto newZone =
                    ConstantFieldOnZone( zone.getFiniteElementDescriptor(), toCopy );
                const auto resu = vect_resu[pos];
                fieldOut->setValueOnZone( newZone, resu );
            }
        }
    };

    void allocateFields( const ThermalLoadPtr< ConstantFieldOnCellsType > &load ) {
        auto typeLoadStd = load->getType();
        typeLoadStd->updateValuePointer();

        _type->allocate( 1 );
        ( *_type )[0] = ( *typeLoadStd )[0];

        _modelName->allocate( 1 );
        ( *_modelName )[0] = _model->getName();

        transferConstantFieldOnCells( load->getImposedField(), _cimpo );
        transferConstantFieldOnCells( load->getMultiplicativeField(), _cmult );
    };

  protected:
    /** @brief Modele */
    ModelPtr _model;
    /** @brief Vecteur Jeveux '.LIGRE' */
    ParallelFiniteElementDescriptorPtr _FEDesc;
    /** @brief Carte '.CIMPO' */
    ConstantFieldOnCellsTypePtr _cimpo;
    /** @brief Carte '.CMULT' */
    ConstantFieldOnCellsRealPtr _cmult;
    /** @brief Vecteur Jeveux '.TYPE' */
    JeveuxVectorChar8 _type;
    /** @brief Vecteur Jeveux '.MODEL.NOMO' */
    JeveuxVectorChar8 _modelName;
    SyntaxSaverPtr _syntax;
    /** @brief node and cell support groups of ParallelThermalLoad */
    VectorString _grpNo, _grpMa;

  public:
    /**
     * @brief Constructeur
     */
    ParallelThermalLoad( void ) = delete;

    /**
     *
     * @brief Constructeur
     */
    ParallelThermalLoad( const ThermalLoadPtr< ConstantFieldOnCellsType > &load,
                         const ModelPtr &model )
        : ParallelThermalLoad( ResultNaming::getNewResultName(), load, model ) {};

    /** @brief Constructor */
    ParallelThermalLoad( const ParallelThermalLoadPtr load, const ModelPtr &model )
        : DataStructure( ResultNaming::getNewResultName(), 8, "CHAR_THER" ),
          _type( getName() + ".TYPE" ),
          _modelName( getName() + ".CHTH.MODEL.NOMO" ),
          _model( model ) {
        _syntax = load->_syntax;
        _grpNo = load->_grpNo;
        _grpMa = load->_grpMa;
        auto pMesh = std::static_pointer_cast< ParallelMesh >( model->getMesh() );
        ConnectionMeshPtr connectionMesh( new ConnectionMesh( pMesh, _grpNo, _grpMa ) );
        ModelPtr connectionModel( new Model( connectionMesh ) );
        connectionModel->setFrom( model );
        ThermalLoadPtr< ConstantFieldOnCellsType > partialMechanicalLoad(
            new ThermalLoad< ConstantFieldOnCellsType >( connectionModel ) );
        partialMechanicalLoad->buildFromSyntax( _syntax );

        _FEDesc = std::make_shared< ParallelFiniteElementDescriptor >(
            getName() + ".CHTH.LIGRE", partialMechanicalLoad->getFiniteElementDescriptor(),
            partialMechanicalLoad->getModel()->getConnectionMesh(), model );
        _cimpo = std::make_shared< ConstantFieldOnCellsType >( getName() + ".CHTH.CIMPO", _FEDesc );
        _cmult = std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CHTH.CMULT", _FEDesc );

        allocateFields( partialMechanicalLoad );
    };

    /**
     * @brief Constructeur
     */
    ParallelThermalLoad( const std::string &name,
                         const ThermalLoadPtr< ConstantFieldOnCellsType > &load,
                         const ModelPtr &model )
        : DataStructure( name, 8, "CHAR_THER" ),
          _FEDesc( std::make_shared< ParallelFiniteElementDescriptor >(
              getName() + ".CHTH.LIGRE", load->getFiniteElementDescriptor(),
              load->getModel()->getConnectionMesh(), model ) ),
          _cimpo(
              std::make_shared< ConstantFieldOnCellsType >( getName() + ".CHTH.CIMPO", _FEDesc ) ),
          _cmult(
              std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CHTH.CMULT", _FEDesc ) ),
          //
          _type( getName() + ".TYPE" ),
          _modelName( getName() + ".CHTH.MODEL.NOMO" ),
          _model( model ) {
        auto typeLoadStd = load->getType();
        typeLoadStd->updateValuePointer();

        _type->allocate( 1 );
        ( *_type )[0] = ( *typeLoadStd )[0];

        _modelName->allocate( 1 );
        ( *_modelName )[0] = model->getName();

        transferConstantFieldOnCells( load->getImposedField(), _cimpo );
        transferConstantFieldOnCells( load->getMultiplicativeField(), _cmult );
    };

    /**
     * @brief Get the finite element descriptor
     */
    ParallelFiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _FEDesc; };

    /**
     * @brief Get the model
     */
    ModelPtr getModel() const { return _model; };

    ConstantFieldOnCellsRealPtr getMultiplicativeField() const { return _cmult; };

    ConstantFieldOnCellsTypePtr getImposedField() const { return _cimpo; };

    /**
     * @brief Function to set rebuild parameters (used for balancing)
     * @param syntax aster syntax
     * @param grpNo cell groups
     * @param syntax node groups
     */
    void setRebuildParameters( SyntaxSaverPtr syntax, const VectorString &grpNo,
                               const VectorString &grpMa ) {
        _syntax = syntax;
        _grpNo = grpNo;
        _grpMa = grpMa;
    };
};

/**
 * @typedef ParallelThermalLoadPtr
 * @brief Pointeur intelligent vers un ParallelThermalLoad
 */
using ParallelThermalLoadReal = ParallelThermalLoad< ConstantFieldOnCellsReal >;

using ParallelThermalLoadFunction = ParallelThermalLoad< ConstantFieldOnCellsChar24 >;

using ParallelThermalLoadRealPtr = std::shared_ptr< ParallelThermalLoadReal >;

using ParallelThermalLoadFunctionPtr = std::shared_ptr< ParallelThermalLoadFunction >;

/** @typedef std::list de ParallelThermalLoad */
using ListParaTherLoadReal = std::list< ParallelThermalLoadRealPtr >;
/** @typedef Iterateur sur une std::list de ParallelThermalLoad */
using ListParaTherLoadRealIter = ListParaTherLoadReal::iterator;
/** @typedef Iterateur constant sur une std::list de ParallelThermalLoad */
using ListParaTherLoadRealCIter = ListParaTherLoadReal::const_iterator;

/** @typedef std::list de ParallelThermalLoad */
using ListParaTherLoadFunction = std::list< ParallelThermalLoadFunctionPtr >;
/** @typedef Iterateur sur une std::list de ParallelThermalLoad */
using ListParaTherLoadFunctionIter = ListParaTherLoadFunction::iterator;
/** @typedef Iterateur constant sur une std::list de ParallelThermalLoad */
using ListParaTherLoadFunctionCIter = ListParaTherLoadFunction::const_iterator;

#endif /* PARALLELTHERMALLOAD_H_ */

#endif /* ASTER_HAVE_MPI */
