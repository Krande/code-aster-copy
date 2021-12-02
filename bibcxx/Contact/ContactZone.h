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

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"
#include "astercxx.h"

class ContactZone : public DataStructure {
  private:
    /** @brief Modele */
    ModelPtr _model;

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
    ContactZone( const std::string name, const ModelPtr model )
        : DataStructure( name, 8, "CHAR_CONT_ZONE" ), _model( model ){};

    /**
     * @brief Constructeur
     */
    ContactZone( const ModelPtr model ) : ContactZone( ResultNaming::getNewResultName(), model ){};

    ModelPtr getModel() const { return _model; }

    BaseMeshPtr getMesh() const { return _model->getMesh(); }
};

/**
 * @typedef ContactZonePtr
 * @brief Pointeur intelligent vers un ContactZone
 */
typedef boost::shared_ptr< ContactZone > ContactZonePtr;

#endif /* CONTACT_ZONE_H_ */
