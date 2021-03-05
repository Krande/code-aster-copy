
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

#include "DataStructures/DataStructure.h"
#include "DataFields/ConstantFieldOnCells.h"
#include "Modeling/Model.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/ParallelFiniteElementDescriptor.h"
#include "Loads/MechanicalLoad.h"
#include "Supervis/ResultNaming.h"

/**
 * @class ParallelMechanicalLoadClass
 * @brief Classe definissant une charge dualisée parallèle
 * @author Nicolas Sellenet
 */
class ParallelMechanicalLoadClass: public DataStructure
{
private:
    void transferConstantFieldOnCells( const ConstantFieldOnCellsRealPtr& fieldIn,
                                ConstantFieldOnCellsRealPtr& fieldOut )
        ;

protected:
    /** @brief Modele */
    ModelPtr                           _model;
    /** @brief Vecteur Jeveux '.LIGRE' */
    ParallelFiniteElementDescriptorPtr _FEDesc;
    /** @brief Carte '.CIMPO' */
    ConstantFieldOnCellsRealPtr        _cimpo;
    /** @brief Carte '.CMULT' */
    ConstantFieldOnCellsRealPtr        _cmult;
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
    ParallelMechanicalLoadClass( const GenericMechanicalLoadPtr& load,
                                    const ModelPtr& model ):
        ParallelMechanicalLoadClass( ResultNaming::getNewResultName(),
                                        load, model )
    {};

    /**
     * @brief Constructeur
     */
    ParallelMechanicalLoadClass( const std::string& name,
                                    const GenericMechanicalLoadPtr& load,
                                    const ModelPtr& model );

    /**
     * @typedef ParallelMechanicalLoadPtr
     * @brief Pointeur intelligent vers un ParallelMechanicalLoad
     */
    typedef boost::shared_ptr< ParallelMechanicalLoadClass > ParallelMechanicalLoadPtr;


    /**
     * @brief Get the finite element descriptor
     */
    ParallelFiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _FEDesc; };

    /**
     * @brief Get the model
     */
    const ModelPtr &getModel() const {
        if ( ( !_model ) || _model->isEmpty() )
            throw std::runtime_error( "Model of current load is empty" );
        return _model;
    };
};

/**
 * @typedef ParallelMechanicalLoadPtr
 * @brief Pointeur intelligent vers un ParallelMechanicalLoadClass
 */
typedef boost::shared_ptr< ParallelMechanicalLoadClass > ParallelMechanicalLoadPtr;

/** @typedef std::list de ParallelMechanicalLoad */
typedef std::list< ParallelMechanicalLoadPtr > ListParaMecaLoad;
/** @typedef Iterateur sur une std::list de ParallelMechanicalLoad */
typedef ListParaMecaLoad::iterator ListParaMecaLoadIter;
/** @typedef Iterateur constant sur une std::list de ParallelMechanicalLoad */
typedef ListParaMecaLoad::const_iterator ListParaMecaLoadCIter;

#endif /* PARALLELMECHANICALLOAD_H_ */

#endif /* ASTER_HAVE_MPI */
