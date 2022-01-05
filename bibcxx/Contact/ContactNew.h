#ifndef CONTACTNEW_H_
#define CONTACTNEW_H_

/**
 * @file ContactNew.h
 * @brief Fichier entete de la class ContactNew
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

#include "Contact/ContactZone.h"
#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"
#include "astercxx.h"

class ContactNew : public DataStructure {
  private:
    /** @brief Modele */
    ModelPtr _model;
    /** @brief Ligel ".CHME.LIGRE" */
    FiniteElementDescriptorPtr _FEDesc;
    /** @brief List of contact zone */
    std::vector< ContactZonePtr > _zones;
    /** @brief Level of verbosity */
    ASTERINTEGER _verbosity;
    /** @brief Friction */
    bool _friction;
    /** @brief Smoothing for normals */
    bool _smoothing;

  public:
    /**
     * @typedef ContactNewPtr
     * @brief Pointeur intelligent vers un ContactNew
     */
    typedef boost::shared_ptr< ContactNew > ContactNewPtr;
    /**
     * @brief Constructeur
     */
    ContactNew() = delete;

    /**
     * @brief Constructeur
     */
    ContactNew( const std::string name, const ModelPtr model );

    /**
     * @brief Constructeur
     */
    ContactNew( const ModelPtr model ) : ContactNew( ResultNaming::getNewResultName(), model ){};

    ModelPtr getModel() const { return _model; }

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

    void hasFriction( const bool &friction ) { _friction = friction; }

    bool hasFriction() const { return _friction; }

    void hasSmoothing( const bool &smoothing ) { _smoothing = smoothing; }

    bool hasSmoothing() const { return _smoothing; }
};

/**
 * @typedef ContactNewPtr
 * @brief Pointeur intelligent vers un ContactNew
 */
typedef boost::shared_ptr< ContactNew > ContactNewPtr;

#endif /* CONTACTNEW_H_ */
