#ifndef CONTACT_ZONE_H_
#define CONTACT_ZONE_H_

/**
 * @file ContactZone.h
 * @brief Fichier entete de la class ContactZone
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

#include "Contact/ContactEnum.h"
#include "Contact/ContactParameters.h"
#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"
#include "astercxx.h"

class ContactZone : public DataStructure {
  private:
    /** @brief Modele */
    ModelPtr _model;
    /** @brief Level of verbosity */
    ASTERINTEGER _verbosity;
    /** @brief Parameter for contact only */
    ContactParameterPtr _contParam;
    /** @brief Parameter for friction only */
    FrictionParameterPtr _fricParam;
    /** @brief Parameter for pairing only */
    PairingParameterPtr _pairParam;
    /** @brief Slave side */
    std::string _slave;
    /** @brief Master side */
    std::string _master;

  public:
    /**
     * @typedef ContactZonePt
     * @brief Pointeur intelligent vers un ContactZone
     */
    typedef boost::shared_ptr< ContactZone > ContactZonePtr;
    /**
     * @brief Constructeur
     */
    ContactZone() = delete;

    /**
     * @brief Constructeur
     */
    ContactZone( const std::string name, const ModelPtr model );

    /**
     * @brief Constructeur
     */
    ContactZone( const ModelPtr model ) : ContactZone( ResultNaming::getNewResultName(), model ){};

    ModelPtr getModel() const { return _model; }

    BaseMeshPtr getMesh() const { return _model->getMesh(); }

    void setVerbosity( const ASTERINTEGER &level ) { _verbosity = level; }

    ASTERINTEGER getVerbosity() const { return _verbosity; }

    bool build();

    ContactParameterPtr getContactParameter() const { return _contParam; };

    FrictionParameterPtr getFrictionParameter() const { return _fricParam; };

    PairingParameterPtr getPairingParameter() const { return _pairParam; };

    void setContactParameter( const ContactParameterPtr contParam ) { _contParam = contParam; };

    void setFrictionParameter( const FrictionParameterPtr fricParam ) { _fricParam = fricParam; };

    void setPairingParameter( const PairingParameterPtr pairParam ) { _pairParam = pairParam; };

    void setSlaveGroupOfCells( const std::string& slave) {_slave = slave;};

    std::string getSlaveGroupOfCells( ) const { return _slave;};

    void setMasterGroupOfCells( const std::string& master) {_master = master;};

    std::string getMasterGroupOfCells( ) const { return _master;};

};

/**
 * @typedef ContactZonePtr
 * @brief Pointeur intelligent vers un ContactZone
 */
typedef boost::shared_ptr< ContactZone > ContactZonePtr;

#endif /* CONTACT_ZONE_H_ */
