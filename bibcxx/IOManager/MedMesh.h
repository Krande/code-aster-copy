#ifndef MEDMESH_H_
#define MEDMESH_H_

/**
 * @file MedMesh.h
 * @brief Fichier entete de la classe MedMesh
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

// aslint: disable=C3010

#include "IOManager/MedCalculationSequence.h"
#include "IOManager/MedFamily.h"
#include "IOManager/MedFilePointer.h"
#include "IOManager/MedJoint.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#ifdef ASTER_HAVE_MED
#include "med.h"
/**
 * @class MedMesh
 * @brief Med mesh interface
 * @author Nicolas Sellenet
 */
class MedMesh {
  private:
    /** @brief mesh name */
    const std::string _name = "";
    /** @brief dimension */
    med_int _dim;
    /** @brief step number */
    med_int _nbstep;
    /** @brief all calcultation sequences */
    std::vector< MedCalculationSequence > _sequences;
    /** @brief all families */
    std::vector< MedFamilyPtr > _families;
    /** @brief med file id */
    const MedFilePointer &_filePtr;
    /** @brief all families */
    std::vector< MedJointPtr > _joints;

  public:
    /**
     * @typedef MedMeshPtr
     * @brief Pointeur intelligent vers un MedMesh
     */
    typedef std::shared_ptr< MedMesh > MedMeshPtr;

    /**
     * @brief Constructor
     * @param filePtr med pointer on file
     * @param name mesh name
     * @param dim dimension
     */
    MedMesh( const MedFilePointer &filePtr, const std::string &name, med_int dim, med_int nbstep );

    /**
     * @brief add a calculation sequence
     * @param numdt step id
     * @param numit iteration id
     * @param dt time step value
     */
    void addSequence( int numdt, int numit, float dt ) {
        if ( _sequences.size() >= _nbstep ) {
            throw std::runtime_error( "Too many sequences defined for mesh " + _name );
        }
        _sequences.emplace_back( numdt, numit, dt );
    };

    /** @brief get cell family from sequence and profile iterator */
    std::vector< med_int > getCellFamilyAtSequence( int numdt, int numit, int iterator ) const;

    /** @brief get cell family from sequence and geotype */
    std::vector< med_int >
    getCellFamilyForGeometricTypeAtSequence( int numdt, int numit,
                                             med_geometry_type geotype ) const;

    /** @brief get cell connectivity from sequence and iterator on field geotype */
    std::vector< med_int > getConnectivityAtSequence( int numdt, int numit, int iterator ) const;

    /** @brief get cell connectivity from sequence and geotype */
    std::vector< med_int >
    getConnectivityForGeometricTypeAtSequence( int numdt, int numit,
                                               med_geometry_type geotype ) const;

    /** @brief get cell number from sequence and iterator on field geotype */
    med_int getCellNumberAtSequence( int numdt, int numit, int iterator ) const;

    /**
     * @brief get all cell number from sequence
     * @return sum of all cell for all geotype in sequence */
    med_int getAllCellNumberAtSequence( int numdt, int numit ) const;

    /** @brief get cell number from sequence and geotype */
    med_int getCellNumberForGeometricTypeAtSequence( int numdt, int numit,
                                                     med_geometry_type geotype ) const;

    /** @brief get geotype from sequence and iterator on field geotype */
    med_geometry_type getCellTypeAtSequence( int numdt, int numit, int iterator ) const;

    /** @brief get type number from sequence */
    med_int getCellTypeNumberAtSequence( int numdt, int numit ) const;

    /** @brief get mesh dimension */
    med_int getDimension() const { return _dim; };

    /** @brief get mesh dimension */
    const std::vector< MedJointPtr > &getJoints() const { return _joints; };

    /** @brief get all families */
    std::vector< MedFamilyPtr > getFamilies() const { return _families; };

    /** @brief get type vector from sequence */
    std::vector< med_int > getGeometricTypesAtSequence( int numdt, int numit ) const;

    /** @brief get global node numbering */
    std::vector< med_int > getGlobalNodeNumberingAtSequence( int numdt, int numit ) const;

    /** @brief get mesh name in med file */
    std::string getName() const { return _name; };

    /** @brief get node family from sequence */
    std::vector< med_int > getNodeFamilyAtSequence( int numdt, int numit ) const;

    /** @brief get node number from sequence */
    med_int getNodeNumberAtSequence( int numdt, int numit ) const;

    /** @brief get node number for geo type */
    med_int getNodeNumberForGeometricType( med_int geoType ) const;

    /** @brief get sequence number for mesh (usualy 1 for aster) */
    int getSequenceNumber() const { return _sequences.size(); };

    /** @brief get sequence (numdt, numit) from id */
    std::vector< med_int > getSequence( int index ) const {
        const auto &curSeq = _sequences[index].getNumDtNumIt();
        return { curSeq.first, curSeq.second };
    };

    /** @brief get split cell number and node id start number from sequence for current proc */
    std::pair< med_int, med_int >
    getSplitCellNumberForGeometricTypeAtSequence( int numdt, int numit,
                                                  med_geometry_type geotype ) const;

    /** @brief get split node number and node id start number from sequence for current proc */
    std::pair< med_int, med_int > getSplitNodeNumberAtSequence( int numdt, int numit ) const;

    /** @brief get coordinates from sequence */
    std::vector< double > readCoordinates( int numdt, int numit ) const;
};

/**
 * @typedef MedMeshPtr
 * @brief Pointeur intelligent vers un MedMesh
 */
typedef std::shared_ptr< MedMesh > MedMeshPtr;

#endif
#endif /* MEDMESH_H_ */
