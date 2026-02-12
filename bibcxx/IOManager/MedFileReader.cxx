/**
 * @file MedFileReader.cxx
 * @brief Implementation de MedFileReader
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

// aslint: disable=C3012

#include "IOManager/MedFileReader.h"

#include "IOManager/MedUtilities.h"
#include "Utilities/Tools.h"

#ifdef ASTER_HAVE_MED
MedFileReader::~MedFileReader() {
    if ( _filePtr.isOpen() )
        _filePtr.close();
};

MedFieldPtr MedFileReader::getField( int index ) const {
    if ( index >= _fields.size() ) {
        throw std::runtime_error( "Field not available" );
    }
    return _fields[index];
};

MedMeshPtr MedFileReader::getMesh( int index ) const {
    if ( index >= _meshes.size() ) {
        throw std::runtime_error( "Mesh not available" );
    }
    return _meshes[index];
};

MedMeshPtr MedFileReader::createMesh( const std::string &name, const int &dim,
                                      const std::string &desc ) {
    for ( const auto &mesh : _meshes ) {
        if ( name == mesh->getName() ) {
            throw std::runtime_error( "Mesh " + name + " already exists in med file" );
        }
    }
    if ( !_filePtr.isOpen() ) {
        throw std::runtime_error( "MED file not open" );
    }
    if ( name.size() > MED_NAME_SIZE ) {
        throw std::runtime_error(
            "Mesh name too long (size limite: " + std::to_string( MED_NAME_SIZE ) + ")" );
    }
    if ( desc.size() > MED_COMMENT_SIZE ) {
        throw std::runtime_error(
            "Mesh description too long (size limite: " + std::to_string( MED_COMMENT_SIZE ) + ")" );
    }
    std::string dtunit( "Sans unite" );
    std::string axisname( "X               Y               Z               " );
    std::string axisunit( "INCONNU         INCONNU         INCONNU         " );
    auto cret = MEDmeshCr( _filePtr.getFileId(), name.c_str(), dim, dim, MED_UNSTRUCTURED_MESH,
                           desc.c_str(), dtunit.c_str(), MED_SORT_ITDT, MED_CARTESIAN,
                           axisname.c_str(), axisunit.c_str() );
    auto ref = _meshes.emplace_back( new MedMesh( _filePtr, name, dim ) );
    return ref;
};

MedFieldPtr MedFileReader::createField( const std::string &name, const MedMeshPtr &mesh,
                                        const VectorString &cmps ) {
    for ( const auto &field : _fields ) {
        if ( name == field->getName() ) {
            throw std::runtime_error( "Field " + name + " already exists in med file" );
        }
    }
    if ( !_filePtr.isOpen() ) {
        throw std::runtime_error( "MED file not open" );
    }
    if ( name.size() > MED_NAME_SIZE ) {
        throw std::runtime_error(
            "Field name too long (size limite: " + std::to_string( MED_NAME_SIZE ) + ")" );
    }
    VectorString cmpUnit;
    for ( int i = 0; i < cmps.size(); ++i ) {
        cmpUnit.push_back( "Sans unite" );
    }
    char *cmpChar = stringVectorToChar( cmps, MED_SNAME_SIZE );
    char *unitChar = stringVectorToChar( cmpUnit, MED_SNAME_SIZE );
    std::string dtunit( "Sans unite" );
    auto cret = MEDfieldCr( _filePtr.getFileId(), name.c_str(), MED_DOUBLE, cmps.size(), cmpChar,
                            unitChar, dtunit.c_str(), mesh->getName().c_str() );
    free( cmpChar );
    free( unitChar );
    auto ref = _fields.emplace_back( new MedField( _filePtr, name, cmps, mesh, 0, {} ) );
    return ref;
};

int MedFileReader::close() { return _filePtr.close(); };

MedFieldPtr MedFileReader::getField( const std::string &name ) const {
    const auto index = _mapFieldNameInt.at( name );
    return _fields[index];
};

std::vector< std::string > MedFileReader::getFieldNames() const {
    VectorString toReturn;
    for ( const auto &keyValue : _mapFieldNameInt )
        toReturn.push_back( keyValue.first );
    return toReturn;
};

int MedFileReader::getFieldNumber() const { return _fields.size(); };

int MedFileReader::getMeshNumber() const { return _meshes.size(); };

int MedFileReader::getProfileNumber() const { return _profiles.size(); };

int MedFileReader::open( const std::filesystem::path &filename, const MedFileAccessType &openType,
                         std::array< int, 3 > version ) {
    _filePtr.open( filename, openType, version );
    readFile();

    return 0;
};

int MedFileReader::openParallel( const std::filesystem::path &filename,
                                 const MedFileAccessType &openType ) {
    _filePtr.openParallel( filename, openType );
    readFile();

    return 0;
};

int MedFileReader::readFile() {
    _fields = std::vector< MedFieldPtr >();
    _profiles = std::vector< MedProfilePtr >();

    const auto fileId = _filePtr.getFileId();

    const auto nbMeshes = MEDnMesh( fileId );
    char meshname[MED_NAME_SIZE + 1] = "";
    char *axisname, *axisunit;
    int count = 0;
    // Read all meshes
    for ( int i = 1; i <= nbMeshes; ++i ) {
        const auto nbAxis = MEDmeshnAxis( fileId, i );
        med_int spacedim = 0, meshdim = 0, nstep = 0;
        med_mesh_type meshtype;
        char description[MED_COMMENT_SIZE + 1] = "";
        char dtunit[MED_SNAME_SIZE + 1] = "";
        med_axis_type axistype;
        med_sorting_type sortingtype;
        axisname = (char *)malloc( ( nbAxis * MED_SNAME_SIZE + 1 ) * sizeof( char ) );
        axisunit = (char *)malloc( ( nbAxis * MED_SNAME_SIZE + 1 ) * sizeof( char ) );
        MEDmeshInfo( fileId, i, meshname, &spacedim, &meshdim, &meshtype, description, dtunit,
                     &sortingtype, &nstep, &axistype, axisname, axisunit );
        if ( meshtype != MED_UNSTRUCTURED_MESH ) {
            std::cout << "Mesh type must be unstructured. Mesh name: " << meshname << std::endl;
            continue;
        }
        auto ref = _meshes.emplace_back( new MedMesh( _filePtr, meshname, spacedim ) );
        ref->reserveSequences( nstep );
        ref->readFromFile();
        _mapMeshNameRank[strip( meshname )] = count;
        free( axisname );
        free( axisunit );
        ++count;
    }

    const auto nbProfile = MEDnProfile( fileId );
    char profilename[MED_NAME_SIZE + 1] = "";
    count = 0;
    // Read all profiles
    for ( int i = 1; i <= nbProfile; ++i ) {
        med_int profilesize = 0;
        MEDprofileInfo( fileId, i, profilename, &profilesize );
        std::string name( profilename );
        _profiles.emplace_back( new MedProfile( name, profilesize ) );
        _mapProfileNameRank[name] = count;
        ++count;
    }

    const auto nbFields = MEDnField( fileId );
    char *componentname, *componentunit;
    char fieldname[MED_NAME_SIZE + 1] = "";
    med_bool localmesh;
    med_field_type fieldtype;
    char dtunit[MED_SNAME_SIZE + 1] = "";
    med_int ncstp = 0;
    count = 0;
    // Read all fields
    for ( int i = 1; i <= nbFields; ++i ) {
        const auto nbCmp = MEDfieldnComponent( fileId, i );
        componentname = (char *)malloc( ( nbCmp * MED_SNAME_SIZE + 1 ) * sizeof( char ) );
        componentunit = (char *)malloc( ( nbCmp * MED_SNAME_SIZE + 1 ) * sizeof( char ) );
        MEDfieldInfo( fileId, i, fieldname, meshname, &localmesh, &fieldtype, componentname,
                      componentunit, dtunit, &ncstp );
        if ( fieldtype != MED_DOUBLE && fieldtype != MED_FLOAT64 && fieldtype != MED_FLOAT32 ) {
            std::cout << "Field type must be float or double. Field name: " << fieldname
                      << std::endl;
            continue;
        }
        const auto cnames = splitChar( componentname, nbCmp, MED_SNAME_SIZE );
        auto curMesh = _meshes[_mapMeshNameRank[strip( meshname )]];
        auto &ref = _fields.emplace_back(
            new MedField( _filePtr, fieldname, cnames, curMesh, ncstp, _profiles ) );

        free( componentname );
        free( componentunit );
        _mapFieldNameInt[strip( fieldname )] = count;
        ++count;
    }
    return 0;
};

#endif
