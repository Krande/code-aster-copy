/**
 * @file CouplingZonePairing.h
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

#include "Contact/ContactParameter.h"
#include "DataFields/FieldOnCells.h"
#include "DataStructures/DataStructure.h"
#include "Meshes/BaseMesh.h"
#include "Meshes/MeshEnum.h"
#include "Meshes/MeshPairing.h"

enum class CouplingMethod {
    Undefined,
    Nitsche,
    Penalization,
    Lagrangian,
};

class CouplingZonePairing : public DataStructure {
    /** Datastructure for pairing */
  protected:
    /** @brief Definition of pairing of two surfaces */
    MeshPairingPtr _meshPairing;

    CouplingMethod _algo;

    ASTERDOUBLE _coef_pena;

    PairingParameterPtr _params;

  public:
    /** @brief No default constructor */
    CouplingZonePairing() = delete;

    /** @brief Constructor with given name */
    CouplingZonePairing( const BaseMeshPtr mesh, const ASTERINTEGER verbosity = 1 );

    /** @brief Get mesh */
    BaseMeshPtr getMesh() const { return _meshPairing->getMesh(); };

    void setSlaveGroupsOfCells( const VectorString &groupName );

    void setMasterGroupsOfCells( const VectorString &groupName );

    /** @brief Compute pairing quantities of all zones */
    ASTERBOOL compute();

    void check( const ModelPtr model ) const;

    ASTERINTEGER getNumberOfPairs() const { return _meshPairing->getNumberOfPairs(); };

    VectorPairLong getListOfPairs() const { return _meshPairing->getListOfPairs(); }

    MapLong getSlaveCellsSurfToVolu() const { return _meshPairing->getSlaveCellsSurfToVolu(); }

    ASTERINTEGER getSlaveCellSurfToVolu( const ASTERINTEGER index ) const {
        return _meshPairing->getSlaveCellSurfToVolu( index );
    }

    VectorOfVectorsReal getIntersectionPoints( const CoordinatesSpace coorSpace ) const;

    VectorOfVectorsReal getIntersectionPoints( const ASTERINTEGER indexPair,
                                               const CoordinatesSpace coorSpace ) const;

    void setMethod( const CouplingMethod algo ) { _algo = algo; }

    CouplingMethod getMethod() const { return _algo; };

    void setCoefficient( const ASTERDOUBLE coef ) { _coef_pena = coef; }

    ASTERDOUBLE getCoefficient() const { return _coef_pena; };

    void setVerbosity( const ASTERINTEGER verbosity );

    void setPairingParameters( const PairingParameterPtr params );

    bool build();
};

using CouplingZonePairingPtr = std::shared_ptr< CouplingZonePairing >;
