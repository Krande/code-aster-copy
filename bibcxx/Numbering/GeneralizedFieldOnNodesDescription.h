#ifndef GENERALIZEDFIELDONNODESDESCRIPTION_H_
#define GENERALIZEDFIELDONNODESDESCRIPTION_H_

/**
 * @file GeneralizedFieldOnNodesDescription.h
 * @brief Fichier entete de la classe GeneralizedFieldOnNodesDescription
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "astercxx.h"

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Supervis/ResultNaming.h"
#include "MemoryManager/JeveuxCollection.h"

/**
 * @class GeneralizedFieldOnNodesDescriptionClass
 * @brief This class describes the structure of dof stored in a field on nodes
 * @author Nicolas Sellenet
 */
class GeneralizedFieldOnNodesDescriptionClass : public DataStructure {
    /** @brief Objet Jeveux '.DESC' */
    JeveuxVectorLong _desc;
    /** @brief Objet Jeveux '.NEQU' */
    JeveuxVectorLong _nequ;
    /** @brief Objet Jeveux '.REFN' */
    JeveuxVectorChar24 _refn;
    /** @brief Objet Jeveux '.DEEQ' */
    JeveuxVectorLong _deeq;
    /** @brief Objet Jeveux '.DELG' */
    JeveuxVectorLong _delg;
    /** @brief Objet Jeveux '.LILI' */
    JeveuxVectorChar24 _lili;
    /** @brief Objet Jeveux '.NUEQ' */
    JeveuxVectorLong _nueq;
    /** @brief Objet Jeveux '.PRNO' */
    JeveuxCollectionLong _prno;
    /** @brief Objet Jeveux '.ORIG' */
    JeveuxCollectionLong _orig;

  public:
    /**
     * @brief Constructeur
     */
    GeneralizedFieldOnNodesDescriptionClass( const JeveuxMemory memType = Permanent )
        : DataStructure( ResultNaming::getNewResultName(), 19, "PROF_GENE", memType ),
          _desc( JeveuxVectorLong( getName() + ".DESC" ) ),
          _nequ( JeveuxVectorLong( getName() + ".NEQU" ) ),
          _refn( JeveuxVectorChar24( getName() + ".REFN" ) ),
          _deeq( JeveuxVectorLong( getName() + ".DEEQ" ) ),
          _delg( JeveuxVectorLong( getName() + ".DELG" ) ),
          _lili( JeveuxVectorChar24( getName() + ".LILI" ) ),
          _nueq( JeveuxVectorLong( getName() + ".NUEQ" ) ),
          _prno( JeveuxCollectionLong( getName() + ".PRNO" ) ),
          _orig( JeveuxCollectionLong( getName() + ".ORIG" ) ){};

    /**
     * @brief Destructor
     */
    ~GeneralizedFieldOnNodesDescriptionClass(){};

    /**
     * @brief Constructeur
     * @param name nom souhaité de la sd (utile pour le GeneralizedFieldOnNodesDescriptionClass
     * d'une sd_resu)
     */
    GeneralizedFieldOnNodesDescriptionClass( const std::string name,
                                                const JeveuxMemory memType = Permanent )
        : DataStructure( name, 19, "PROF_GENE", memType ),
          _desc( JeveuxVectorLong( getName() + ".DESC" ) ),
          _nequ( JeveuxVectorLong( getName() + ".NEQU" ) ),
          _refn( JeveuxVectorChar24( getName() + ".REFN" ) ),
          _deeq( JeveuxVectorLong( getName() + ".DEEQ" ) ),
          _delg( JeveuxVectorLong( getName() + ".DELG" ) ),
          _lili( JeveuxVectorChar24( getName() + ".LILI" ) ),
          _nueq( JeveuxVectorLong( getName() + ".NUEQ" ) ),
          _prno( JeveuxCollectionLong( getName() + ".PRNO" ) ),
          _orig( JeveuxCollectionLong( getName() + ".ORIG" ) ){};
};
typedef boost::shared_ptr< GeneralizedFieldOnNodesDescriptionClass >
    GeneralizedFieldOnNodesDescriptionPtr;

#endif /* GENERALIZEDFIELDONNODESDESCRIPTION_H_ */
