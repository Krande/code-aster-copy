#ifndef PARALLELMESH_H_
#define PARALLELMESH_H_

/**
 * @file ParallelMesh.h
 * @brief Fichier entete de la classe ParallelMesh
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

#include <set>

#ifdef ASTER_HAVE_MPI

#include "MemoryManager/NamesMap.h"
#include "Meshes/BaseMesh.h"
#include "Supervis/ResultNaming.h"

/**
 * @class ParallelMesh
 * @brief Cette classe decrit un maillage Aster parallèle
 * @author Nicolas Sellenet
 */
class ParallelMesh : public BaseMesh {
  private:
    typedef std::set< std::string > SetOfString;
    typedef SetOfString::iterator SetOfStringIter;
    typedef SetOfString::const_iterator SetOfStringCIter;

    /** @brief All groups of nodes (parallel mesh) */
    NamesMapChar24 _globalGroupOfNodes;
    /** @brief Set of all groups of nodes (parallel mesh) */
    SetOfString _setOfAllGON;
    /** @brief All groups of cells (parallel mesh) */
    NamesMapChar24 _globalGroupOfCells;
    /** @brief Set of all groups of cells (parallel mesh) */
    SetOfString _setOfAllGOE;
    /** @brief Identify outer nodes */
    JeveuxVectorLong _outerNodes;
    /** @brief Identify outer cells */
    JeveuxVectorLong _outerCells;
    /** @brief Global numbering */
    JeveuxVectorLong _globalNumbering;
    /** @brief List of opposite domain */
    JeveuxVectorLong _listOfOppositeDomain;
    /** @brief List of joints */
    std::map< ASTERINTEGER, std::pair< JeveuxVectorLong, JeveuxVectorLong > > _joints;
    /** @brief Global to local node number */
    std::unordered_map< ASTERINTEGER, ASTERINTEGER > _global2localMap;

    void _buildGlobal2LocalMap() {
        getLocalToGlobalMapping()->updateValuePointer();
        ASTERINTEGER nloc = getLocalToGlobalMapping()->size();

        _global2localMap.reserve( nloc );
        for ( auto j = 0; j < nloc; j++ )
            _global2localMap[( *getLocalToGlobalMapping() )[j]] = j;
    };

  public:
    /**
     * @typedef ParallelMeshPtr
     * @brief Pointeur intelligent vers un ParallelMesh
     */
    typedef std::shared_ptr< ParallelMesh > ParallelMeshPtr;

    /**
     * @brief Constructeur
     */
    ParallelMesh() : ParallelMesh( ResultNaming::getNewResultName() ){};

    /**
     * @brief Constructeur
     */
    ParallelMesh( const std::string &name )
        : BaseMesh( name, "MAILLAGE_P" ),
          _globalGroupOfNodes( getName() + ".PAR_GRPNOE" ),
          _globalGroupOfCells( getName() + ".PAR_GRPMAI" ),
          _outerNodes( getName() + ".NOEX" ),
          _outerCells( getName() + ".MAEX" ),
          _globalNumbering( getName() + ".NULOGL" ),
          _listOfOppositeDomain( getName() + ".DOMJOINTS" ){};

    /**
     * @brief Get the JeveuxVector for outer subdomain nodes
     * @return _outerNodes
     */
    const JeveuxVectorLong getNodesRank() const { return _outerNodes; };

    /**
     * @brief Get the JeveuxVector for outer subdomain cells
     * @return _outerCells
     */
    const JeveuxVectorLong getCellsRank() const { return _outerCells; };

    bool hasGroupOfCells( const std::string &name, const bool local = false ) const;

    bool hasGroupOfNodes( const std::string &name, const bool local = false ) const;

    VectorString getGroupsOfCells( const bool local = false ) const;

    VectorString getGroupsOfNodes( const bool local = false ) const;

    VectorLong getCells( const std::string name = "" ) const;

    /**
     * @brief Return list of nodes
     * @param name name of group (if empty all the nodes)
     * @param localNumbering node id in local or global numbering
     * @param same_rank keep or not the nodes owned by the current domain
     * @return list of Nodes
     */

    VectorLong getNodes( const std::string name = std::string(), const bool localNumbering = true,
                         const ASTERINTEGER same_rank = PythonBool::None ) const;

    /**
     * @brief Get inner nodes
     * @return list of node ids
     */
    VectorLong getInnerNodes() const {
        return this->getNodes( std::string(), true, PythonBool::True );
    };

    /**
     * @brief Get outer nodes
     * @return list of node ids
     */
    VectorLong getOuterNodes() const {
        return this->getNodes( std::string(), true, PythonBool::False );
    };

    /**
     * @brief Get inner nodes
     * @return list of node ids
     */
    VectorLong getInnerCells() const;

    /**
     * @brief Get outer nodes
     * @return list of node ids
     */
    VectorLong getOuterCells() const;

    /**
     * @brief Get the mapping between local and global numbering of nodes
     * @return JeveuxVector of the indirection
     */
    const JeveuxVectorLong getLocalToGlobalMapping() const {return _globalNumbering;};

    /**
     * @brief Returns the nodes indexes of a group of cells
     * @param name name of group of cells
     * @param local node id in local or global numbering
     * @param same_rank keep or not the nodes owned by the current domain
     * @return list of Nodes
     */
    VectorLong getNodesFromCells( const std::string name, const bool localNumbering = true,
                                  const ASTERINTEGER same_rank = PythonBool::None ) const;

    VectorLong getNodesFromCells( const VectorLong &cells, const bool localNumbering = true,
                                  const ASTERINTEGER same_rank = PythonBool::None ) const;

    /**
     * @brief Fonction permettant de savoir si un maillage est parallel
     * @return retourne true si le maillage est parallel
     */
    virtual bool isParallel() const { return true; };

    /**
     * @brief Read a MED ParallelMesh file (already partitioned mesh)
     * @return retourne true si tout est ok
     */
    bool readPartitionedMedFile( const std::string &fileName );

    /**
     * @brief Get the local number of a node from its local number
     * @return Return the local number if the node if present on the subdomain ; 
     * otherwise raise an exception
     */
    const ASTERINTEGER globalToLocalNodeId(const ASTERINTEGER);

    bool updateGlobalGroupOfNodes( void );

    bool updateGlobalGroupOfCells( void );

    bool build();
};

/**
 * @typedef ParallelMeshPtr
 * @brief Pointeur intelligent vers un ParallelMesh
 */
typedef std::shared_ptr< ParallelMesh > ParallelMeshPtr;

#endif /* ASTER_HAVE_MPI */

#endif /* PARALLELMESH_H_ */
