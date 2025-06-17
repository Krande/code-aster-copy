/**
 * @file ContactNew.h
 * @brief Header of class ContactNew
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

#pragma once

#include "astercxx.h"

#include "Contact/ContactZone.h"
#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"

class ContactNew : public DSWithCppPickling {
  protected:
    /** @brief Model */
    ModelPtr _model;

    /** @brief Ligrel ".CONT.LIGRE" (virtual cells) */
    FiniteElementDescriptorPtr _FEDesc;

    /** @brief List of contact zones */
    std::vector< ContactZonePtr > _zones;

    /** @brief Level of verbosity */
    ASTERINTEGER _verbosity;

  protected:
    /** @brief Main constructor */
    ContactNew( const std::string name, const ModelPtr model, const std::string type );

  private:
    /** @brief Get global dimension of space */
    ASTERINTEGER getSpaceDime() const;

  public:
    using ContactNewPtr = std::shared_ptr< ContactNew >;

    /** @brief No default constructor */
    ContactNew() = delete;

    /** @brief Constructor with given name */
    ContactNew( const std::string name, const ModelPtr model )
        : ContactNew( name, model, "CHAR_CONT" ) {};

    /** @brief Constructor with automatic name */
    ContactNew( const ModelPtr model ) : ContactNew( ResultNaming::getNewResultName(), model ) {};

    /** @brief restricted constructor (Set) and method (Get) to support pickling */
    ContactNew( const py::tuple &tup );
    py::tuple _getState() const;

    /** @brief Get mesh */
    BaseMeshPtr getMesh() const { return _model->getMesh(); }

    /** @brief Get model */
    ModelPtr getModel() const { return _model; }

    /** @brief Get Finite Element Descriptor */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _FEDesc; }

    /** @brief Set verbosity */
    void setVerbosity( const ASTERINTEGER &level );

    /** @brief Get verbosity */
    ASTERINTEGER getVerbosity() const { return _verbosity; }

    /** @brief Append contact zone */
    void appendContactZone( const ContactZonePtr zone );

    /** @brief Get number of contact zones */
    ASTERINTEGER getNumberOfContactZones() const { return _zones.size(); }

    /** @brief Get contact zone */
    ContactZonePtr getContactZone( const ASTERINTEGER &zone_id ) const {
        return _zones.at( zone_id );
    }

    /** @brief Get pair geometry of zone */
    MeshPairingPtr getMeshPairing( const ASTERINTEGER &zone_id ) {
        return _zones.at( zone_id )->getMeshPairing();
    };

    /** @brief Clear pairing result of zone */
    void clearPairing( const ASTERINTEGER &zone_id ) {
        _zones.at( zone_id )->getMeshPairing()->clearResult();
    };

    /** @brief Set coordinates */
    void setCoordinates( const MeshCoordinatesFieldPtr &coor ) {
        for ( auto &zone : _zones ) {
            zone->setCoordinates( coor );
        }
    };

    /** @brief Get all contact zones */
    std::vector< ContactZonePtr > getContactZones() const { return _zones; }

    /** @brief Get all slaves nodes */
    VectorLong getSlaveNodes() const;

    /** @brief Get all slaves cells */
    VectorLong getSlaveCells() const;

    /** @brief Set/unset friction flag everywhere */
    void enableFriction( const bool &friction );

    /** @brief Detect if friction on one of the contact zone */
    bool hasFriction() const;

    /** @brief Set/unset smoothing of normals flag everywhere */
    void enableSmoothing( const bool &smoothing );

    /** @brief Detect if smoothing of normals on one of the contact zone */
    bool hasSmoothing() const;

    /** @brief Get number of intersection points on all zones */
    VectorLong getNumberOfIntersectionPoints() const;

    /** @brief Get number of intersection points on a contact zone */
    VectorLong getNumberOfIntersectionPoints( const ASTERINTEGER &indexZone ) const;

    /** @brief Builder from Fortran part */
    virtual bool build();

    /** @brief To know if ContactNew (or child) is Parallel */
    bool isParallel() const { return false; };
};

using ContactNewPtr = std::shared_ptr< ContactNew >;

class FrictionNew : public ContactNew {

  public:
    using FrictionNewPtr = std::shared_ptr< FrictionNew >;

    /** @brief No default constructor */
    FrictionNew() = delete;

    /** @brief Constructor with given name */
    FrictionNew( const std::string name, const ModelPtr model )
        : ContactNew( name, model, "CHAR_FROT" ) {};

    /** @brief Constructor with automatic name */
    FrictionNew( const ModelPtr model ) : FrictionNew( ResultNaming::getNewResultName(), model ) {};

    /** @brief restricted constructor (Set) and method (Get) to support pickling */
    FrictionNew( const py::tuple &tup );
    py::tuple _getState() const;

    /** @brief Builder from Fortran part */
    bool build() {
        AS_ASSERT( hasFriction() );
        return ContactNew::build();
    };
};

using FrictionNewPtr = std::shared_ptr< FrictionNew >;
