#ifndef FIELDONNODESDESCRIPTION_H_
#define FIELDONNODESDESCRIPTION_H_

/**
 * @file FieldOnNodesDescription.h
 * @brief Fichier entete de la classe FieldOnNodesDescription
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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
#include "DataStructures/DataStructureNaming.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "MemoryManager/NamesMap.h"
#include "Supervis/ResultNaming.h"

/**
 * @class FieldOnNodesDescription
 * @brief This class describes the structure of dof stored in a field on nodes
 * @author Nicolas Sellenet
 */
class FieldOnNodesDescription : public DataStructure {
    /** @brief Objet Jeveux '.PRNO' */
    JeveuxCollectionLong _componentsOnNodes;
    /** @brief Objet Jeveux '.LILI' */
    NamesMapChar24 _namesOfGroupOfCells;
    /** @brief Objet Jeveux '.NUEQ' */
    JeveuxVectorLong _indexationVector;
    /** @brief Objet Jeveux '.DEEQ' */
    JeveuxVectorLong _nodeAndComponentsNumberFromDOF;

  public:
    /**
     * @brief Constructeur
     * @param name nom souhaité de la sd (utile pour le FieldOnNodesDescription d'une
     * sd_resu)
     */
    FieldOnNodesDescription( const std::string name );

    /**
     * @brief Constructeur
     */
    FieldOnNodesDescription() : FieldOnNodesDescription( DataStructureNaming::getNewName() ){};

    /**
     * @brief Destructor
     */
    ~FieldOnNodesDescription(){};

    /**
     * @brief Surcharge de l'operateur =
     */
    bool operator==( FieldOnNodesDescription &toCompare );

    /**
     * @brief Returns a vector of information of the numbering
     */
    const JeveuxVectorLong getNodeAndComponentsNumberFromDOF() const {
        return _nodeAndComponentsNumberFromDOF;
    }

    /**
     * @brief Returns a vector with node index for each DOFs
     */
    VectorLong getNodesFromDOF() const;

    /**
     * @brief Returns number of DOFs
     */
    ASTERINTEGER getNumberOfDofs() const;

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
     */
    bool updateValuePointers();
};

typedef boost::shared_ptr< FieldOnNodesDescription > FieldOnNodesDescriptionPtr;

#endif /* FIELDONNODESDESCRIPTION_H_ */
