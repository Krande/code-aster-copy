#ifndef INCOMPLETEMESH_H_
#define INCOMPLETEMESH_H_

/**
 * @file IncompleteMesh.h
 * @brief Fichier entete de la classe IncompleteMesh
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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
#include "Meshes/Mesh.h"
#include "Supervis/ResultNaming.h"

/**
 * @class IncompleteMesh
 * @brief Class of an incomplete mesh read in parallel from medcoupling
 * @author Nicolas Sellenet
 */
class IncompleteMesh : public Mesh {

    VectorLong _nodeRange;
    VectorOfVectorsLong _cellRange;
    VectorLong _nodeFamily;
    VectorLong _cellFamily;
    std::vector< VectorString > _nodeFamGroups;
    std::vector< VectorString > _cellFamGroups;

  public:
    /**
     * @typedef IncompleteMeshPtr
     * @brief Pointeur intelligent vers un IncompleteMesh
     */
    typedef std::shared_ptr< IncompleteMesh > IncompleteMeshPtr;

    /**
     * @brief Constructeur
     */
    IncompleteMesh() : IncompleteMesh( ResultNaming::getNewResultName() ) {};

    /**
     * @brief Constructeur
     */
    IncompleteMesh( const std::string &name ) : Mesh( name, "MAILLAGE_I" ) {};

    void addFamily( int id, VectorString groups );

    bool debugCheckFromBaseMesh( BaseMeshPtr & ) const;

    ASTERINTEGER getDimension() const;

    VectorLong getCellFamily() const { return _cellFamily; };

    std::vector< VectorString > getCellFamilyGroups() const { return _cellFamGroups; };

    const VectorOfVectorsLong &getCellRange() const { return _cellRange; };

    std::vector< VectorString > getNodeFamilyGroups() const { return _nodeFamGroups; };

    VectorLong getNodeFamily() const { return _nodeFamily; };

    VectorLong getNodesFromGroup( const std::string &grpName ) const;

    const VectorLong &getNodeRange() const { return _nodeRange; };

    bool isIncomplete() const { return true; };

    bool isParallel() const { return false; };

    void setCellFamily( const VectorLong &cf );

    void setCellRange( const VectorOfVectorsLong &range ) { _cellRange = range; };

    void setNodeFamily( const VectorLong &nf ) { _nodeFamily = nf; };

    void setNodeRange( const VectorLong &range ) { _nodeRange = range; };
};

/**
 * @typedef IncompleteMeshPtr
 * @brief Pointeur intelligent vers un IncompleteMesh
 */
typedef std::shared_ptr< IncompleteMesh > IncompleteMeshPtr;

#endif /* ASTER_HAVE_MPI */

#endif /* INCOMPLETEMESH_H_ */
