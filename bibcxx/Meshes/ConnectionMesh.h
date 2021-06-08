#ifndef CONNECTIONMESH_H_
#define CONNECTIONMESH_H_

/**
 * @file ConnectionMesh.h
 * @brief Fichier entete de la classe
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

#include "astercxx.h"

#ifdef ASTER_HAVE_MPI

#include "Meshes/BaseMesh.h"
#include "Meshes/ParallelMesh.h"
#include "Supervis/ResultNaming.h"

/**
 * @class ConnectionMeshClass
 * @brief Cette classe decrit un maillage partiel reconstruit a partir d'une liste de groupe de
 * noeuds
 * @author Nicolas Sellenet
 */
class ConnectionMeshClass : public BaseMeshClass {
  private:
    typedef JeveuxCollection< ASTERINTEGER, NamesMapChar24 > JeveuxCollectionLongNamePtr;
    /** @brief Base ParallelMesh */
    ParallelMeshPtr _pMesh;
    /** @brief id of node in local numbering */
    JeveuxVectorLong _localNumbering;
    /** @brief id of node in global numbering */
    JeveuxVectorLong _globalNumbering;
    /** @brief number of owner proc for each nodes */
    JeveuxVectorLong _owner;
    /** @brief id of cell in local numbering */
    JeveuxVectorLong _cellsLocalNumbering;
    /** @brief number of owner proc for each cells */
    JeveuxVectorLong _cellsOwner;


    VectorLong getCellsGlobalNumbering( const JeveuxVectorLong& rankOfCells ) const;

  public:
    /**
     * @typedef ConnectionMeshPtr
     * @brief Pointeur intelligent vers un ConnectionMeshClass
     */
    typedef boost::shared_ptr< ConnectionMeshClass > ConnectionMeshPtr;

    /**
     * @brief Constructeur
     */
    ConnectionMeshClass( const ParallelMeshPtr &mesh, const VectorString &groupsOfNodes,
        const VectorString &groupsOfCells )
        : ConnectionMeshClass( ResultNaming::getNewResultName(), mesh,
                            groupsOfNodes, groupsOfCells ){};

    /**
     * @brief Constructeur
     */
    ConnectionMeshClass( const std::string &name, const ParallelMeshPtr &mesh,
        const VectorString &groupsOfNodes, const VectorString &groupsOfCells );

    const JeveuxVectorLong &getGlobalNumbering() const { return _globalNumbering; };

    const JeveuxVectorLong &getLocalNumbering() const { return _localNumbering; };

    const JeveuxVectorLong &getOwner() const { return _owner; };

    const JeveuxVectorLong &getCellsLocalNumbering() const { return _cellsLocalNumbering; };

    const JeveuxVectorLong &getCellsOwner() const { return _cellsOwner; };

    const ParallelMeshPtr &getParallelMesh() const { return _pMesh; };

    VectorString getGroupsOfCells( ) const;

    VectorString getGroupsOfNodes( ) const;

    bool hasGroupOfCells( const std::string &name) const;

    bool hasGroupOfNodes( const std::string &name) const;

    const VectorLong getCells( const std::string name ) const;

    const VectorLong getCells(  ) const
    {
        return getCells( std::string() );
    };

    /**
     * @brief Fonction permettant de savoir si un maillage est partiel
     * @return retourne true si le maillage est partiel
     */
    virtual bool isConnection() const { return true; };
};

/**
 * @typedef ConnectionMeshPtr
 * @brief Pointeur intelligent vers un ConnectionMeshClass
 */
typedef boost::shared_ptr< ConnectionMeshClass > ConnectionMeshPtr;

#endif /* ASTER_HAVE_MPI */

#endif /* CONNECTIONMESH_H_ */
