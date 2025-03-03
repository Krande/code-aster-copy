#ifndef MESHREADER_H_
#define MESHREADER_H_

/**
 * @file MeshReader.h
 * @brief Fichier entete de la classe MeshReader
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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

#include "Meshes/Mesh.h"

#include <filesystem>

/**
 * @class MeshReader
 * @brief Med file interface
 * @author Nicolas Sellenet
 */
class MeshReader {
  private:
  public:
    /**
     * @typedef MeshReaderPtr
     * @brief Pointeur intelligent vers un MeshReader
     */
    typedef std::shared_ptr< MeshReader > MeshReaderPtr;

    /** @brief Constructor */
    MeshReader() {};

    ~MeshReader() {};
#ifdef ASTER_HAVE_MED

    MeshPtr readFromMedFile( const std::filesystem::path &filename );
    // void testPerf();
#endif
};

/**
 * @typedef MeshReaderPtr
 * @brief Pointeur intelligent vers un MeshReader
 */
typedef std::shared_ptr< MeshReader > MeshReaderPtr;

#endif /* MESHREADER_H_ */
