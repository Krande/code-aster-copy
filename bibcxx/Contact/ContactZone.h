#ifndef CONTACT_ZONE_H_
#define CONTACT_ZONE_H_

/**
 * @file ContactZone.h
 * @brief Header of class ContactZone
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

#include "astercxx.h"

#include "Contact/ContactEnum.h"
#include "Contact/ContactParameter.h"
#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/MeshPairing.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"

class ContactZone : public DSWithCppPickling {
  private:
    /** @brief Model */
    ModelPtr _model;
    /** @brief Level of verbosity */
    ASTERINTEGER _verbosity;
    /** @brief Parameters for contact */
    ContactParameterPtr _contParam;
    /** @brief Parameters for friction  */
    FrictionParameterPtr _fricParam;
    /** @brief Parameters for pairing */
    PairingParameterPtr _pairParam;
    /** @brief Definition of pairing of two surfaces */
    MeshPairingPtr _meshPairing;
    /** @brief Check direction of normal */
    bool _checkNormal;
    /** @brief  Smoothing of normal */
    bool _smoothing;

  public:
    using ContactZonePtr = std::shared_ptr< ContactZone >;

    /** @brief Constructor with given name */
    ContactZone( const std::string name );

    /** @brief Constructor with automatic name */
    ContactZone() : ContactZone( ResultNaming::getNewResultName() ) {};

    /** @brief restricted constructor (Set) and method (Get) to support pickling */
    ContactZone( const py::tuple &tup );
    py::tuple _getState() const;

    /** @brief Get model */
    ModelPtr getModel() const { return _model; }

    /** @brief Get mesh */
    BaseMeshPtr getMesh() const { return _model->getMesh(); }

    /** @brief Set verbosity */
    void setVerbosity( const ASTERINTEGER &level );

    /** @brief Get verbosity */
    ASTERINTEGER getVerbosity() const { return _verbosity; }

    /** @brief Get parameters for contact */
    ContactParameterPtr getContactParameter() const { return _contParam; };

    /** @brief Get parameters for friction */
    FrictionParameterPtr getFrictionParameter() const { return _fricParam; };

    /** @brief Get parameters for pairing */
    PairingParameterPtr getPairingParameter() const { return _pairParam; };

    /** @brief Set parameters for contact */
    void setContactParameter( const ContactParameterPtr contParam ) { _contParam = contParam; };

    /** @brief Set parameters for friction */
    void setFrictionParameter( const FrictionParameterPtr fricParam ) { _fricParam = fricParam; };

    /** @brief Set parameters for pairing */
    void setPairingParameter( const PairingParameterPtr pairParam ) { _pairParam = pairParam; };

    /** @brief Set group of slave cells */
    void setSlaveGroupOfCells( const std::string &groupName ) {
        _meshPairing->setSlaveGroupOfCells( groupName );
    }

    /** @brief Set group of master cells */
    void setMasterGroupOfCells( const std::string &groupName ) {
        _meshPairing->setMasterGroupOfCells( groupName );
    }

    /** @brief Get master cells */
    const VectorLong &getMasterCells() const { return _meshPairing->getMasterCells(); };
    VectorLong &getMasterCells() {
        return const_cast< VectorLong & >( std::as_const( *this ).getMasterCells() );
    }

    /** @brief Get slave cells */
    VectorLong getSlaveCells() const { return _meshPairing->getSlaveCells(); }

    /** @brief Set excluded groups of slave cells */
    void setExcludedSlaveGroupOfCells( const VectorString &groupsName ) {
        _meshPairing->setExcludedSlaveGroupOfCells( groupsName );
    }

    /** @brief Set excluded groups of slave nodes */
    void setExcludedSlaveGroupOfNodes( const VectorString &groupsName ) {
        _meshPairing->setExcludedSlaveGroupOfNodes( groupsName );
    }

    /** @brief Set/get check normals */
    void checkNormals( const bool &checkNormal ) { _checkNormal = checkNormal; }
    bool checkNormals() const { return _checkNormal; }

    /** @brief Get master nodes */
    VectorLong getMasterNodes() const { return _meshPairing->getMasterNodes(); };

    /** @brief Get slave nodes */
    VectorLong getSlaveNodes() const { return _meshPairing->getSlaveNodes(); };

    /** @brief Get volume slave cells linked to all surfacic slave cells */
    MapLong getSlaveCellsSurfToVolu() const { return _meshPairing->getSlaveCellsSurfToVolu(); };

    /** @brief Get volume slave cell linked to a surfacic slave cells */
    ASTERINTEGER getSlaveCellSurfToVolu( const ASTERINTEGER &cellIndex ) const {
        return _meshPairing->getSlaveCellSurfToVolu( cellIndex );
    };

    /** @brief Set coordinates */
    void setCoordinates( const MeshCoordinatesFieldPtr &coor ) {
        _meshPairing->setCoordinates( coor );
    };

    /** @brief Set/unset friction for this zone */
    void enableFriction( const bool &friction ) { _fricParam->enableFriction( friction ); };

    /** @brief Detect if friction for this zone */
    bool hasFriction() const { return _fricParam->hasFriction(); };

    /** @brief Set/unset normal smoothing for this zone */
    void enableSmoothing( const bool &smoothing ) { _smoothing = smoothing; };

    /** @brief Detect if normal smoothing for this zone */
    bool hasSmoothing() const { return _smoothing; };

    /** @brief Compute pairing of zone */
    bool pairing( ASTERDOUBLE &dist_pairing, ASTERDOUBLE &pair_tole );

    /** @brief Get pairing of surface meshes */
    MeshPairingPtr getMeshPairing() { return _meshPairing; };

    /** @brief Builder from Fortran part */
    bool build( const ModelPtr model );

    const std::string &getSlaveGroupOfCells() const {
        return _meshPairing->getSlaveGroupOfCells();
    };

    const std::string &getMasterGroupOfCells() const {
        return _meshPairing->getMasterGroupOfCells();
    };
};

using ContactZonePtr = std::shared_ptr< ContactZone >;

#endif /* CONTACT_ZONE_H_ */
