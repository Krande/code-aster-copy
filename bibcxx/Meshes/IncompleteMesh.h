#ifndef INCOMPLETEMESH_H_
#define INCOMPLETEMESH_H_

/**
 * @file IncompleteMesh.h
 * @brief Fichier entete de la classe IncompleteMesh
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

#include <set>

#ifdef ASTER_HAVE_MPI

#include "MemoryManager/NamesMap.h"
#include "Meshes/Mesh.h"
#include "Supervis/ResultNaming.h"

/**
 * @class IncompleteMesh
 * @brief Cette classe decrit un maillage Aster parallèle
 * @author Nicolas Sellenet
 */
class IncompleteMesh : public Mesh {
  public:
    /**
     * @typedef IncompleteMeshPtr
     * @brief Pointeur intelligent vers un IncompleteMesh
     */
    typedef std::shared_ptr< IncompleteMesh > IncompleteMeshPtr;

    /**
     * @brief Constructeur
     */
    IncompleteMesh() : IncompleteMesh( ResultNaming::getNewResultName() ){};

    /**
     * @brief Constructeur
     */
    IncompleteMesh( const std::string &name ) : Mesh( name, "MAILLAGE_I" ){};

    bool isIncomplete() const { return true; };

    bool isParallel() const { return false; };
};

/**
 * @typedef IncompleteMeshPtr
 * @brief Pointeur intelligent vers un IncompleteMesh
 */
typedef std::shared_ptr< IncompleteMesh > IncompleteMeshPtr;

#endif /* ASTER_HAVE_MPI */

#endif /* INCOMPLETEMESH_H_ */
