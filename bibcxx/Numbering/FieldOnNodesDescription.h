#ifndef FIELDONNODESDESCRIPTION_H_
#define FIELDONNODESDESCRIPTION_H_

/**
 * @file FieldOnNodesDescription.h
 * @brief Fichier entete de la classe FieldOnNodesDescription
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

#include "DataStructures/DataStructure.h"
#include "DataStructures/DataStructureNaming.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "MemoryManager/NamesMap.h"
#include "Meshes/BaseMesh.h"
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
    /** @brief Mesh (only in c++) */
    BaseMeshPtr _mesh;

  public:
    /** @typedef FieldOnNodesDescriptionPtr */
    typedef std::shared_ptr< FieldOnNodesDescription > FieldOnNodesDescriptionPtr;

    /**
     * @brief Constructeur
     * @param name nom souhaité de la sd (utile pour le FieldOnNodesDescription d'une
     * sd_resu)
     */
    FieldOnNodesDescription( const std::string name, const BaseMeshPtr mesh = nullptr,
                             const std::string type = "PROF_CHNO" );

    /**
     * @brief Constructeur
     */
    FieldOnNodesDescription( const BaseMeshPtr mesh = nullptr )
        : FieldOnNodesDescription( DataStructureNaming::getNewName(), mesh ){};

    /**
     * @brief Destructor
     */
    ~FieldOnNodesDescription(){};

    /**
     * @brief Surcharge de l'operateur =
     */
    bool operator==( FieldOnNodesDescription &toCompare );

    bool operator!=( FieldOnNodesDescription &toCompare ) { return !( *this == toCompare ); }

    /**
     * @brief Returns a vector with node index for each DOFs
     */
    VectorLong getNodesFromDOF() const;

    /**
     * @brief Returns number of DOFs
     */
    ASTERINTEGER getNumberOfDofs() const;

    /**
     * @brief Return list of DOFs
     * @param sameRank True: Use only owned nodes / False: Use all nodes
     * @param list_cmp empty: Use all cmp / keep only cmp given
     * @param groupsOfCells empty: Use all nodes / keep only nodes given
     */
    VectorLong getDOFs( const bool sameRank = false, const VectorString &list_cmp = {},
                        const VectorLong &list_nodes = {} ) const;

    /**
     * @brief Returns a vector with node index and component name for each DOFs
     */
    VectorPairLong getNodesAndComponentsNumberFromDOF( const bool local = true ) const;

    PairLong getNodeAndComponentNumberFromDOF( const ASTERINTEGER dof,
                                               const bool local = true ) const;

    /**
     * @brief Returns a vector with node index and component name for each DOFs
     */
    std::vector< std::pair< ASTERINTEGER, std::string > >
    getNodesAndComponentsFromDOF( const bool local = true ) const;
    std::pair< ASTERINTEGER, std::string >
    getNodeAndComponentFromDOF( const ASTERINTEGER dof, const bool local = true ) const;

    /**
     * @brief Maps between node id and name of components to DOF
     */
    std::map< PairLong, ASTERINTEGER >
    getDOFsFromNodesAndComponentsNumber( const bool local = true ) const;

    std::map< std::pair< ASTERINTEGER, std::string >, ASTERINTEGER >
    getDOFsFromNodesAndComponents( const bool local = true ) const;

    /**
     * @brief Get componants
     */
    SetString getComponents() const;
    SetLong getComponentsNumber() const;

    /**
     * @brief Maps between name of components and the number
     */
    std::map< std::string, ASTERINTEGER > getComponentsName2Number() const;
    std::map< ASTERINTEGER, std::string > getComponentsNumber2Name() const;

    void setMesh( const BaseMeshPtr mesh );

    BaseMeshPtr getMesh() const { return _mesh; };

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
     */
    void updateValuePointers();
};

typedef std::shared_ptr< FieldOnNodesDescription > FieldOnNodesDescriptionPtr;

#endif /* FIELDONNODESDESCRIPTION_H_ */
