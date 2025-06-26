#ifndef MEDFIELD_H_
#define MEDFIELD_H_

/**
 * @file MedField.h
 * @brief Fichier entete de la classe MedField
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

#include "IOManager/MedCalculationSequence.h"
#include "IOManager/MedFilePointer.h"
#include "IOManager/MedMesh.h"
#include "IOManager/MedProfile.h"
#include "IOManager/MedTypes.h"
#include "IOManager/MedVector.h"

#include <array>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#ifdef ASTER_HAVE_MED
#include "med.h"
// aslint: disable=C3010
// aslint: disable=C3012

/**
 * @class MedField
 * @brief Med field interface
 * @author Nicolas Sellenet
 */
class MedField {
  private:
    /** @brief field name */
    std::string _name = "";
    /** @brief component list */
    std::vector< std::string > _componentName;
    /** @brief support mesh */
    const MedMeshPtr _mesh;
    /** @brief iteration number */
    med_int _nbIt = 0, _nbCmp;
    /** @brief sequence vector */
    std::vector< MedCalculationSequence > _sequences;
    /** @brief map between step id/iteration id and position in _sequences */
    std::map< med_int, std::map< med_int, med_int > > _numdtNumitToSeq;
    /** @brief id of last added sequence */
    int curSeq = 1;
    const MedFilePointer &_filePtr;
    /** @brief all profiles in med file */
    std::vector< MedProfilePtr > _profiles;
    /** @brief map between profile name and profile rank */
    std::map< std::string, int > _mapProfileNameRank;

  public:
    /**
     * @typedef MedFieldPtr
     * @brief Pointeur intelligent vers un MedField
     */
    typedef std::shared_ptr< MedField > MedFieldPtr;

    /**
     * @brief Constructor
     * @param filePtr MedFilePointer of field
     * @param name field name in med file
     * @param componentName vector of component name
     * @param mesh pointer on MedMesh
     * @param nbIt number of iterations
     * @param nbCmp component number
     * @param profiles vector of all profiles in med file
     */
    MedField( const MedFilePointer &filePtr, const std::string &name,
              const std::vector< std::string > &componentName, const MedMeshPtr &mesh, med_int nbIt,
              med_int nbCmp, std::vector< MedProfilePtr > profiles )
        : _name( name ),
          _componentName( componentName ),
          _mesh( mesh ),
          _nbIt( nbIt ),
          _nbCmp( nbCmp ),
          _filePtr( filePtr ),
          _profiles( profiles ) {
        for ( int i = 0; i < _profiles.size(); ++i )
            _mapProfileNameRank[_profiles[i]->getName()] = i;
    };

    /**
     * @brief Add a calculation sequence
     * @param numdt step id
     * @param numit iteration id
     * @param dt time step value
     */
    void addSequence( int numdt, int numit, float dt ) {
        _sequences.emplace_back( numdt, numit, dt );
        _numdtNumitToSeq[numdt][numit] = curSeq;
        ++curSeq;
    };

    /** @brief Get vector of all entity type and geometry type in calculation sequence */
    std::vector< med_int > getAllSupportEntitiesAtSequence( int numdt, int numit ) const;

    /** @brief Get component name */
    std::vector< std::string > getComponentName() const { return _componentName; };

    /** @brief Get component number */
    int getComponentNumber() const { return _nbCmp; };

    /** @brief Get profile number from sequence, med entity and med geotype */
    int getProfileNumberAtSequenceOnEntity( int numdt, int numit, int ent,
                                            med_geometry_type geo ) const;

    /** @brief Get name */
    std::string getName() const { return _name; };

    /** @brief Get sequence from id */
    std::vector< med_int > getSequence( int index ) const {
        const auto &curSeq = _sequences[index].getNumDtNumIt();
        return { curSeq.first, curSeq.second };
    };

    /** @brief Get time from id */
    med_float getTime( int index ) const { return _sequences[index].getDt(); };

    /** @brief Get sequence number */
    int getSequenceNumber() const { return _sequences.size(); };

    /** @brief Get values on cells from sequence, med entity and a list of geotype */
    MedVectorPtr getValuesAtSequenceOnCellTypesList( int numdt, int numit,
                                                     std::vector< med_geometry_type > ) const;

    /** @brief Get values on nodes from sequence */
    MedVectorPtr getValuesAtSequenceOnNodes( int numdt, int numit ) const;

    /** @brief Get values on nodes from sequence, geotype and profile id */
    std::vector< double > getValuesAtSequenceOnEntityAndProfile( int numdt, int numit, int ent,
                                                                 med_geometry_type geo,
                                                                 int profileit ) const;
};

/**
 * @typedef MedFieldPtr
 * @brief Pointeur intelligent vers un MedField
 */
typedef std::shared_ptr< MedField > MedFieldPtr;

#endif
#endif /* MEDFIELD_H_ */
