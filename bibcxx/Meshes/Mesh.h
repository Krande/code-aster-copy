#ifndef MESH_H_
#define MESH_H_

/**
 * @file Mesh.h
 * @brief Fichier entete de la classe Mesh
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

#include "Meshes/BaseMesh.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Supervis/ResultNaming.h"

/**
 * @class Mesh
 * @brief Cette classe decrit un maillage Aster
 */
class Mesh : public BaseMesh {
  protected:
    /**
     * @brief Constructeur
     */
    Mesh( const std::string name, const std::string type ) : BaseMesh( name, type ){};

  public:
    /**
     * @typedef MeshPtr
     * @brief Pointeur intelligent vers un Mesh
     */
    typedef std::shared_ptr< Mesh > MeshPtr;

    /**
     * @brief Constructeur
     */
    Mesh() : BaseMesh( ResultNaming::getNewResultName(), "MAILLAGE" ){};

    /**
     * @brief Constructeur
     */
    Mesh( const std::string name ) : BaseMesh( name, "MAILLAGE" ){};

    bool hasGroupOfCells( const std::string &name, const bool local = false ) const;

    bool hasGroupOfNodes( const std::string &name, const bool local = false ) const;

    VectorString getGroupsOfCells( const bool local = false ) const;

    VectorString getGroupsOfNodes( const bool local = false ) const;

    VectorLong getCells( const std::string name = "" ) const;

    /**
     * @brief Return list of nodes
     * @param name name of group (if empty all the nodes)
     * @param local node id in local or global numbering
     * @param same_rank keep or not the nodes owned by the current domain
     * @return list of Nodes
     */
    VectorLong getNodes( const std::string name = std::string(), const bool localNumbering = true,
                         const ASTERINTEGER same_rank = PythonBool::None ) const;

    /**
     * @brief Returns the nodes indexes of a group of cells
     * @param name name of group of cells
     * @param local node id in local or global numbering
     * @param same_rank keep or not the nodes owned by the current domain
     * @return list of nodes indexes
     */
    VectorLong getNodesFromCells( const std::string, const bool localNumbering = true,
                                  const ASTERINTEGER same_rank = PythonBool::None ) const;

    VectorLong getNodesFromCells( const VectorLong& cells, const bool localNumbering = true,
                                  const ASTERINTEGER same_rank = PythonBool::None ) const;

    /**
     * @brief Get inner nodes
     * @return list of node ids
     */
    VectorLong getInnerNodes() const { return this->getNodes(); };

    /**
     * @brief Get inner nodes
     * @return list of node ids
     */
    const JeveuxVectorLong getNodesRank() const {
        return JeveuxVectorLong( getName() + ".NOEX",
                                 VectorLong( getNumberOfNodes(), getMPIRank() ) );
    };

    /**
     * @brief Get inner cells
     * @return list of cells ids
     */
    const JeveuxVectorLong getCellsRank() const {
        return JeveuxVectorLong( getName() + ".MAEX",
                                 VectorLong( getNumberOfCells(), getMPIRank() ) );
    };

    bool isQuadratic() const;

    /**
     * @brief Read a Aster Mesh file
     * @return retourne true si tout est ok
     */
    bool readAsterFile( const std::string &fileName );

    /**
     * @brief Read a Gibi Mesh file
     * @return retourne true si tout est ok
     */
    bool readGibiFile( const std::string &fileName );

    /**
     * @brief Read a Gmsh Mesh file
     * @return retourne true si tout est ok
     */
    bool readGmshFile( const std::string &fileName );

};

/**
 * @typedef MeshPtr
 * @brief Pointeur intelligent vers un Mesh
 */
typedef std::shared_ptr< Mesh > MeshPtr;

#endif /* MESH_H_ */
