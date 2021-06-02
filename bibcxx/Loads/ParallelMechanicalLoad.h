
#include "astercxx.h"

#ifdef ASTER_HAVE_MPI

#ifndef PARALLELMECHANICALLOAD_H_
#define PARALLELMECHANICALLOAD_H_

/**
 * @file ParallelMechanicalLoad.h
 * @brief Fichier entete de la classe ParallelMechanicalLoad
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

#include "DataFields/ConstantFieldOnCells.h"
#include "DataStructures/DataStructure.h"
#include "DataFields/ConstantFieldOnCells.h"
#include "Modeling/Model.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/ParallelFiniteElementDescriptor.h"
#include "Loads/MechanicalLoad.h"
#include "Supervis/ResultNaming.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Meshes/MeshExplorer.h"


/**
 * @class ParallelMechanicalLoadClass
 * @brief Classe definissant une charge dualisée parallèle
 * @author Nicolas Sellenet
 */
template< typename ConstantFieldOnCellsType>
class ParallelMechanicalLoadClass: public DataStructure
{
  public:
    typedef boost::shared_ptr< ConstantFieldOnCellsType > ConstantFieldOnCellsTypePtr;

private:
    void transferConstantFieldOnCells( const ConstantFieldOnCellsTypePtr& fieldIn,
                                ConstantFieldOnCellsTypePtr& fieldOut )
    {
        const auto& toKeep = _FEDesc->getDelayedElementsToKeep();

        std::string savedName( "" );
        fieldOut->allocate( Permanent, fieldIn );
        const auto sizeFieldIn = (*fieldIn).size();
        const auto vect_resu = fieldIn->getAllValues();

        fieldIn->updateValuePointers();
        for( int pos = 0; pos < sizeFieldIn; ++pos )
        {
            const auto& zone = fieldIn->getZoneDescription( pos );
            const auto& curFEDesc = zone.getFiniteElementDescriptor();
            if( curFEDesc->getName() != savedName && savedName != "" )
            {
                std::string a(
                    "Different FiniteElementDescriptor in one ConstantFieldOnCells is not allowed");
                throw std::runtime_error( a );
            }
            savedName = curFEDesc->getName();

            const auto& listCells = zone.getListOfCells();
            VectorLong toCopy;
            toCopy.reserve(listCells.size());
            for( const auto& num : listCells )
            {
                if( toKeep[ -num - 1 ] != 1 )
                    toCopy.push_back( toKeep[ -num - 1 ] );
            }

            if( toCopy.size() != 0 )
            {
                const auto newZone =
                    ConstantFieldOnZone( zone.getFiniteElementDescriptor(), toCopy );
                const auto resu = vect_resu[pos];
                fieldOut->setValueOnZone( newZone, resu );
            }
        }

    };

protected:
    /** @brief Modele */
    ModelPtr                           _model;
    /** @brief Vecteur Jeveux '.LIGRE' */
    ParallelFiniteElementDescriptorPtr _FEDesc;
    /** @brief Carte '.CIMPO' */
    ConstantFieldOnCellsTypePtr        _cimpo;
    /** @brief Carte '.CMULT' */
    ConstantFieldOnCellsTypePtr        _cmult;
    /** @brief Vecteur Jeveux '.TYPE' */
    JeveuxVectorChar8                  _type;
    /** @brief Vecteur Jeveux '.MODEL.NOMO' */
    JeveuxVectorChar8                  _modelName;

public:
    /**
     * @brief Constructeur
     */
    ParallelMechanicalLoadClass( void ) = delete;

    /**
     *
     * @brief Constructeur
     */
    ParallelMechanicalLoadClass(
        const MechanicalLoadPtr< ConstantFieldOnCellsType >& load,
        const ModelPtr& model ):
        ParallelMechanicalLoadClass( ResultNaming::getNewResultName(), load, model )
    {};

    /**
     * @brief Constructeur
     */
    ParallelMechanicalLoadClass( const std::string& name,
                                const MechanicalLoadPtr< ConstantFieldOnCellsType >& load,
                                const ModelPtr& model ):
    DataStructure( name, 8, "CHAR_MECA" ),
    _FEDesc( boost::make_shared<ParallelFiniteElementDescriptorClass>
                    ( getName() + ".CHME.LIGRE", load->getFiniteElementDescriptor(),
                      load->getModel()->getConnectionMesh(), model ) ),
    _cimpo(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".CHME.CIMPO", _FEDesc )),
    _cmult(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".CHME.CMULT", _FEDesc )),
    _type( getName() + ".TYPE" ),
    _modelName( getName() + ".CHME.MODEL.NOMO" ),
    _model( model )
    {
        auto typeLoadStd = load->getType();
        typeLoadStd->updateValuePointer();

        _type->allocate( Permanent, 1 );
        (*_type)[0] = (*typeLoadStd)[0];
        _modelName->allocate( Permanent, 1 );
        (*_modelName)[0] = model->getName();

        transferConstantFieldOnCells(
            load->getMechanicalLoadDescription()->getImpositionValues(), _cimpo );
        transferConstantFieldOnCells(
            load->getMechanicalLoadDescription()->getCoefficient(), _cmult );
    };

    /**
     * @brief Get the finite element descriptor
     */
    ParallelFiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _FEDesc; };

    /**
     * @brief Get the model
     */
    ModelPtr getModel() const {
        if ( ( !_model ) || _model->isEmpty() )
            throw std::runtime_error( "Model of current load is empty" );
        return _model;
    };

    typedef boost::shared_ptr< ParallelMechanicalLoadClass > ParallelMechanicalLoadPtr;

};

/**
 * @typedef ParallelMechanicalLoadPtr
 * @brief Pointeur intelligent vers un ParallelMechanicalLoadClass
 */
typedef ParallelMechanicalLoadClass< ConstantFieldOnCellsRealClass >
    ParallelMechanicalLoadRealClass;

typedef ParallelMechanicalLoadClass< ConstantFieldOnCellsChar24Class >
    ParallelMechanicalLoadFunctionClass;

typedef boost::shared_ptr< ParallelMechanicalLoadRealClass > ParallelMechanicalLoadRealPtr;

typedef boost::shared_ptr< ParallelMechanicalLoadFunctionClass > ParallelMechanicalLoadFunctionPtr;

/** @typedef std::list de ParallelMechanicalLoad */
typedef std::list< ParallelMechanicalLoadRealPtr > ListParaMecaLoadReal;
/** @typedef Iterateur sur une std::list de ParallelMechanicalLoad */
typedef ListParaMecaLoadReal::iterator ListParaMecaLoadRealIter;
/** @typedef Iterateur constant sur une std::list de ParallelMechanicalLoad */
typedef ListParaMecaLoadReal::const_iterator ListParaMecaLoadRealCIter;

/** @typedef std::list de ParallelMechanicalLoad */
typedef std::list< ParallelMechanicalLoadFunctionPtr > ListParaMecaLoadFunction;
/** @typedef Iterateur sur une std::list de ParallelMechanicalLoad */
typedef ListParaMecaLoadFunction::iterator ListParaMecaLoadFunctionIter;
/** @typedef Iterateur constant sur une std::list de ParallelMechanicalLoad */
typedef ListParaMecaLoadFunction::const_iterator ListParaMecaLoadFunctionCIter;

#endif /* PARALLELMECHANICALLOAD_H_ */

#endif /* ASTER_HAVE_MPI */
