#ifndef MEDJOINT_H_
#define MEDJOINT_H_

/**
 * @file MedJoint.h
 * @brief Fichier entete de la classe MedJoint
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "astercxx.h"

#ifdef ASTER_HAVE_MED
#include "med.h"

#include "IOManager/MedFilePointer.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

/**
 * @class MedJoint
 * @brief Med joint interface
 * @author Nicolas Sellenet
 */
class MedJoint {
  private:
    struct CorrespondenceDescription {
        med_entity_type _localEntityType;
        med_entity_type _remoteEntityType;
        med_geometry_type _localGeoType;
        med_geometry_type _remoteGeoType;
        med_int _entityNumber;

        CorrespondenceDescription( med_entity_type localEntityType, med_geometry_type localGeoType,
                                   med_entity_type remoteEntityType,
                                   med_geometry_type remoteGeoType, med_int entityNumber )
            : _localEntityType( localEntityType ),
              _remoteEntityType( remoteEntityType ),
              _localGeoType( localGeoType ),
              _remoteGeoType( remoteGeoType ),
              _entityNumber( entityNumber ) {};
    };

    /** @brief med file id */
    const MedFilePointer &_filePtr;
    /** @brief joint name */
    std::string _meshName;
    /** @brief joint id in mesh */
    int _id;
    /** @brief joint name */
    std::string _name;
    /** @brief joint description */
    std::string _description;
    /** @brief opposite domain */
    int _domain;
    /** @brief opposite mesh name */
    std::string _remoteMeshName;
    /** @brief step number */
    int _stepNb;
    /** @brief corresp number */
    int _corrNb;
    /** @brief correspondence step description */
    std::vector< std::pair< med_int, med_int > > _stepDesc;
    /** @brief correspondence number by step */
    std::vector< med_int > _correspNbVector;
    /** @brief correspondence description */
    std::vector< std::vector< CorrespondenceDescription > > _corrDesc;

  public:
    /**
     * @typedef MedJointPtr
     * @brief Pointeur intelligent vers un MedJoint
     */
    typedef std::shared_ptr< MedJoint > MedJointPtr;

    /** @brief Constructor */
    MedJoint( const MedFilePointer &filePtr, const std::string &meshName, int id,
              const std::string &name, const std::string &description, int domain,
              const std::string &remoteMeshName, int stepNb, int corrNb );

    /** @brief Get correspondence */
    std::vector< med_int > getCorrespondence( int step, int corresp ) const;

    /** @brief Get corresp number */
    int getCorrespondenceNumber() const { return _corrNb; };

    /** @brief Get family id in med file */
    med_int getId() const { return _id; };

    /** @brief Get family name */
    const std::string &getName() const { return _name; };

    /** @brief Get opposite domain */
    int getOppositeDomain() const { return _domain; };

    /** @brief Get step number */
    int getStepNumber() const { return _stepNb; };
};

/**
 * @typedef MedJointPtr
 * @brief Pointeur intelligent vers un MedJoint
 */
typedef std::shared_ptr< MedJoint > MedJointPtr;

#endif
#endif /* MEDJOINT_H_ */
