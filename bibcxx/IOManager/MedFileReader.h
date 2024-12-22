#ifndef MEDFILEREADER_H_
#define MEDFILEREADER_H_

/**
 * @file MedFileReader.h
 * @brief Fichier entete de la classe MedFileReader
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

// aslint: disable=C3012

#include "IOManager/MedField.h"
#include "IOManager/MedFilePointer.h"
#include "IOManager/MedMesh.h"
#include "IOManager/MedProfile.h"

#include <filesystem>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#ifdef ASTER_HAVE_MED
#include "med.h"
/**
 * @class MedFileReader
 * @brief Med file interface
 * @author Nicolas Sellenet
 */
class MedFileReader {
  private:
    /** @brief med file pointer */
    MedFilePointer _filePtr;
    /** @brief vector of all fields in file */
    std::vector< MedFieldPtr > _fields;
    /** @brief vector of all meshes in file */
    std::vector< MedMeshPtr > _meshes;
    /** @brief vector of all profile in file */
    std::vector< MedProfilePtr > _profiles;
    /** @brief map field name -> index in _fields */
    std::map< std::string, int > _mapFieldNameInt;
    /** @brief map mesh name -> index in _meshes */
    std::map< std::string, int > _mapMeshNameRank;
    /** @brief map profile name -> index in _profiles */
    std::map< std::string, int > _mapProfileNameRank;

    int readFile();

  public:
    /**
     * @typedef MedFileReaderPtr
     * @brief Pointeur intelligent vers un MedFileReader
     */
    typedef std::shared_ptr< MedFileReader > MedFileReaderPtr;

    /** @brief Constructor */
    MedFileReader() {};

    ~MedFileReader();

    /** @brief close file */
    int close();

    /** @brief get field from name */
    MedFieldPtr getField( const std::string &name ) const;

    /** @brief get field from index */
    MedFieldPtr getField( int index ) const { return _fields[index]; };

    /** @brief get all field names */
    std::vector< std::string > getFieldNames() const;

    /** @brief get number of field */
    int getFieldNumber() const;

    /** @brief get number of mesh */
    int getMeshNumber() const;

    /** @brief get mesh from index */
    MedMeshPtr getMesh( int index ) const { return _meshes[index]; };

    /** @brief get number of profile */
    int getProfileNumber() const;

    /** @brief med parallel open of file */
    int openParallel( const std::filesystem::path &filename, const MedFileAccessType &openType );

    /** @brief med open of file */
    int open( const std::filesystem::path &filename, const MedFileAccessType &openType );
};

/**
 * @typedef MedFileReaderPtr
 * @brief Pointeur intelligent vers un MedFileReader
 */
typedef std::shared_ptr< MedFileReader > MedFileReaderPtr;

#endif
#endif /* MEDFILEREADER_H_ */
