/**
 * @file ContactNew.h
 * @brief Fichier entete de la class ContactNew
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

#pragma once

#include "astercxx.h"

#include "Contact/ContactZone.h"
#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"

class ContactNew : public DataStructure {
  protected:
    /** @brief Modele */
    ModelPtr _model;
    /** @brief Ligel ".CHME.LIGRE" */
    FiniteElementDescriptorPtr _FEDesc;
    /** @brief List of contact zone */
    std::vector< ContactZonePtr > _zones;
    /** @brief Level of verbosity */
    ASTERINTEGER _verbosity;

    /**
     * @brief Constructeur
     */
    ContactNew( const std::string name, const ModelPtr model, const std::string type );

  public:
    typedef std::vector< std::pair< ASTERINTEGER, ASTERINTEGER > > VectorLongPairs;

    /**
     * @typedef ContactNewPtr
     * @brief Pointeur intelligent vers un ContactNew
     */
    typedef std::shared_ptr< ContactNew > ContactNewPtr;

    /**
     * @brief Constructeur
     */
    ContactNew() = delete;

    /**
     * @brief Constructeur
     */
    ContactNew( const std::string name, const ModelPtr model )
        : ContactNew( name, model, "CHAR_CONT" ) {};

    /**
     * @brief Constructeur
     */
    ContactNew( const ModelPtr model ) : ContactNew( ResultNaming::getNewResultName(), model ) {};

    /**
     * @brief Get Model
     */
    ModelPtr getModel() const { return _model; }

    /**
     * @brief Get FiniteElementDescriptor
     */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _FEDesc; }

    BaseMeshPtr getMesh() const { return _model->getMesh(); }

    void appendContactZone( const ContactZonePtr zone ) { _zones.push_back( zone ); }

    ASTERINTEGER getNumberOfContactZones() const { return _zones.size(); }

    ContactZonePtr getContactZone( const ASTERINTEGER &zone_id ) const {
        return _zones.at( zone_id );
    }

    std::vector< ContactZonePtr > getContactZones() const { return _zones; }

    void setVerbosity( const ASTERINTEGER &level ) { _verbosity = level; }

    ASTERINTEGER getVerbosity() const { return _verbosity; }

    bool build();

    void enableFriction( const bool &friction );

    bool hasFriction() const;

    void enableSmoothing( const bool &smoothing );

    bool hasSmoothing() const;

    VectorLong getSlaveNodes() const;

    VectorLong getSlaveCells() const;
};

/**
 * @typedef ContactNewPtr
 * @brief Pointeur intelligent vers un ContactNew
 */
using ContactNewPtr = std::shared_ptr< ContactNew >;

class FrictionNew : public ContactNew {

  public:
    /**
     * @typedef FrictionNewPtr
     * @brief Pointeur intelligent vers un FrictionNew
     */
    typedef std::shared_ptr< FrictionNew > FrictionNewPtr;

    /**
     * @brief Constructeur
     */
    FrictionNew() = delete;

    /**
     * @brief Constructeur
     */
    FrictionNew( const std::string name, const ModelPtr model )
        : ContactNew( name, model, "CHAR_FROT" ) {};

    /**
     * @brief Constructeur
     */
    FrictionNew( const ModelPtr model ) : FrictionNew( ResultNaming::getNewResultName(), model ) {};

    bool build() {
        AS_ASSERT( hasFriction() );

        return ContactNew::build();
    };
};

/**
 * @typedef FrictionNewPtr
 * @brief Pointeur intelligent vers un FrictionNew
 */
using FrictionNewPtr = std::shared_ptr< FrictionNew >;
