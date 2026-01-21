/**
 * @file CouplingPairing.h
 * @section LICENCE
 *   Copyright (C) 1991 - 2026  EDF www.code-aster.org
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

#pragma once

#include "DataFields/FieldOnCells.h"
#include "DataStructures/DataStructure.h"
#include "Meshes/MeshEnum.h"
#include "Meshes/MeshPairing.h"
#include "Modeling/Model.h"

enum class CouplingMethod {
    Undefined,
    Nitsche,
    Penalization,
};

class CouplingPairing : public DataStructure {
    /** Datastructure for pairing */
  protected:
    /** @brief Model */
    ModelPtr _model;

    /** @brief Finite element descriptor for virtual elements of coupling */
    FiniteElementDescriptorPtr _fed;

    /** @brief Level of verbosity */
    ASTERINTEGER _verbosity;

    /** @brief Definition of pairing of two surfaces */
    MeshPairingPtr _meshPairing;

    CouplingMethod _algo;

    ASTERDOUBLE _coef_pena;

  protected:
    /** @brief Create virtual elements for coupling */
    void createVirtualElemForCoupling( MapLong &cplElemType,
                                       const JeveuxContiguousCollectionLong meshConnectivity,
                                       std::vector< VectorLong > &listCplElem,
                                       std::vector< VectorPairLong > &listCplType,
                                       SetLong &slaveNodePaired, SetLong &slaveCellPaired );

    /** @brief Get index of coupling cell */
    ASTERINTEGER getCplCellType( const std::string &slavCellTypeName,
                                 const std::string &mastCellTypeName ) const;

    /* Check mesh */
    void check() const;

  public:
    /** @brief No default constructor */
    CouplingPairing() = delete;

    /** @brief Constructor with given name */
    CouplingPairing( const std::string name, const ModelPtr model,
                     const ASTERINTEGER verbosity = 1 );

    /** @brief Constructor with automatic name */
    CouplingPairing( const ModelPtr model, const ASTERINTEGER verbosity = 1 )
        : CouplingPairing( ResultNaming::getNewResultName(), model, verbosity ) {};

    /** @brief Get model */
    ModelPtr getModel() const { return _model; };

    /** @brief Get mesh */
    BaseMeshPtr getMesh() const { return _model->getMesh(); };

    /** @brief Compute pairing quantities of all zones */
    ASTERBOOL compute();

    /** @brief Set group of slave cells */
    void setSlaveGroupOfCells( const std::string &groupName ) {
        _meshPairing->setSlaveGroupOfCells( groupName );
    }

    /** @brief Set group of master cells */
    void setMasterGroupOfCells( const std::string &groupName ) {
        _meshPairing->setMasterGroupOfCells( groupName );
    }

    /** @brief Build Finite Element Descriptor from pairing */
    virtual void buildFiniteElementDescriptor();

    /** @brief Get Finite Element Descriptor from pairing */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _fed; };

    /** @brief Set verbosity */
    void setVerbosity( const ASTERINTEGER &level );

    /** @brief Get verbosity */
    ASTERINTEGER getVerbosity() const { return _verbosity; }

    ASTERINTEGER getNumberOfPairs() const { return _meshPairing->getNumberOfPairs(); };

    FieldOnCellsRealPtr getPairingField() const;

    void setMethod( const CouplingMethod algo ) { _algo = algo; }

    void setCoefficient( const ASTERDOUBLE coef ) { _coef_pena = coef; }
};

using CouplingPairingPtr = std::shared_ptr< CouplingPairing >;
