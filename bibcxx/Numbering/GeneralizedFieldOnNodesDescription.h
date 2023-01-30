#ifndef GENERALIZEDFIELDONNODESDESCRIPTION_H_
#define GENERALIZEDFIELDONNODESDESCRIPTION_H_

/**
 * @file GeneralizedFieldOnNodesDescription.h
 * @brief Fichier entete de la classe GeneralizedFieldOnNodesDescription
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "astercxx.h"

#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "Numbering/FieldOnNodesDescription.h"
#include "Supervis/ResultNaming.h"

/**
 * @class GeneralizedFieldOnNodesDescription
 * @brief This class describes the structure of dof stored in a field on nodes
 * @author Nicolas Sellenet
 */
class GeneralizedFieldOnNodesDescription : public FieldOnNodesDescription {
    /** @brief Objet Jeveux '.DESC' */
    JeveuxVectorLong _desc;
    /** @brief Objet Jeveux '.NEQU' */
    JeveuxVectorLong _nequ;
    /** @brief Objet Jeveux '.REFN' */
    JeveuxVectorChar24 _refn;
    /** @brief Objet Jeveux '.DELG' */
    JeveuxVectorLong _delg;
    /** @brief Objet Jeveux '.ORIG' */
    JeveuxCollectionLong _orig;

  public:
    /**
     * @brief Constructeur
     * @param name nom souhaité de la sd (utile pour le GeneralizedFieldOnNodesDescription
     * d'une sd_resu)
     */
    GeneralizedFieldOnNodesDescription( const std::string name )
        : FieldOnNodesDescription( name, "PROF_GENE" ),
          _desc( JeveuxVectorLong( getName() + ".DESC" ) ),
          _nequ( JeveuxVectorLong( getName() + ".NEQU" ) ),
          _refn( JeveuxVectorChar24( getName() + ".REFN" ) ),
          _delg( JeveuxVectorLong( getName() + ".DELG" ) ),
          _orig( JeveuxCollectionLong( getName() + ".ORIG" ) ){};

    /**
     * @brief Constructeur
     */
    GeneralizedFieldOnNodesDescription()
        : GeneralizedFieldOnNodesDescription( ResultNaming::getNewResultName() ){};

    /**
     * @brief Destructor
     */
    ~GeneralizedFieldOnNodesDescription(){};
};
typedef std::shared_ptr< GeneralizedFieldOnNodesDescription > GeneralizedFieldOnNodesDescriptionPtr;

#endif /* GENERALIZEDFIELDONNODESDESCRIPTION_H_ */
