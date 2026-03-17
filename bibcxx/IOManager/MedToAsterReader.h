#ifndef MEDTOASTERREADER_H_
#define MEDTOASTERREADER_H_

/**
 * @file MedToAsterReader.h
 * @brief Fichier entete de la classe MedToAsterReader
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2026  EDF www.code-aster.org
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

#include "IOManager/MedFileReader.h"
#include "Meshes/IncompleteMesh.h"
#include "Meshes/Mesh.h"
#include "Meshes/ParallelMesh.h"

#include <filesystem>

/**
 * @class MedToAsterReader
 * @brief Med file interface
 * @author Nicolas Sellenet
 */
class MedToAsterReader {
  private:
#ifdef ASTER_HAVE_MED
    void _readMesh( BaseMeshPtr, MedFileReader &, const std::string &, int verbosity = 0 );
#endif

  public:
    /**
     * @typedef MedToAsterReaderPtr
     * @brief Pointeur intelligent vers un MedToAsterReader
     */
    typedef std::shared_ptr< MedToAsterReader > MedToAsterReaderPtr;

    /** @brief Constructor */
    MedToAsterReader() {};

    ~MedToAsterReader() {};

#ifdef ASTER_HAVE_MED
    void readMeshFromMedFile( MeshPtr &, const std::filesystem::path &filename,
                              const std::string &meshName = "", int verbosity = 0 );

#ifdef ASTER_HAVE_MPI
    void readIncompleteMeshFromMedFile( IncompleteMeshPtr &, const std::filesystem::path &filename,
                                        const std::string &meshName = "", int verbosity = 0 );

    void readParallelMeshFromMedFile( ParallelMeshPtr &, const std::filesystem::path &filename,
                                      const std::string &meshName = "", int verbosity = 0 );
#endif
#endif
};

/**
 * @typedef MedToAsterReaderPtr
 * @brief Pointeur intelligent vers un MedToAsterReader
 */
typedef std::shared_ptr< MedToAsterReader > MedToAsterReaderPtr;

#endif /* MEDTOASTERREADER_H_ */
