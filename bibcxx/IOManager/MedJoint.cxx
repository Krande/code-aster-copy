/**
 * @file MedJoint.cxx
 * @brief Implementation de MedJoint
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

#include "IOManager/MedJoint.h"

#ifdef ASTER_HAVE_MED

MedJoint::MedJoint( const MedFilePointer &filePtr, const std::string &meshName, int id,
                    const std::string &name, const std::string &description, int domain,
                    const std::string &remoteMeshName, int stepNb, int corrNb )
    : _filePtr( filePtr ),
      _meshName( meshName ),
      _id( id ),
      _name( name ),
      _description( description ),
      _domain( domain ),
      _remoteMeshName( remoteMeshName ),
      _stepNb( stepNb ),
      _corrNb( corrNb ) {

    for ( int stepId = 1; stepId <= _stepNb; ++stepId ) {
        med_int numdt = -1, numit = -1, ncorrespondence = -1;
        MEDsubdomainComputingStepInfo( _filePtr.getFileId(), meshName.c_str(), _name.c_str(),
                                       stepId, &numdt, &numit, &ncorrespondence );
        _stepDesc.push_back( { numdt, numit } );
        _correspNbVector.push_back( ncorrespondence );

        std::vector< CorrespondenceDescription > toPush;
        for ( int corrId = 1; corrId <= ncorrespondence; ++corrId ) {
            med_entity_type localentitype, remoteentitype;
            med_geometry_type localgeotype, remotegeotype;
            med_int nentitycor = -1;
            MEDsubdomainCorrespondenceSizeInfo(
                _filePtr.getFileId(), meshName.c_str(), _name.c_str(), numdt, numit, corrId,
                &localentitype, &localgeotype, &remoteentitype, &remotegeotype, &nentitycor );

            toPush.push_back( CorrespondenceDescription(
                localentitype, localgeotype, remoteentitype, remotegeotype, nentitycor ) );
        }
        _corrDesc.push_back( toPush );
    }
};

std::vector< med_int > MedJoint::getCorrespondence( int step, int corresp ) const {
    const auto &dtIt = _stepDesc[step - 1];
    const auto &correspDesc = _corrDesc[step - 1][corresp - 1];
    std::vector< med_int > toReturn( 2 * correspDesc._entityNumber, 0 );

    auto ier = MEDsubdomainCorrespondenceRd(
        _filePtr.getFileId(), _meshName.c_str(), _name.c_str(), dtIt.first, dtIt.second,
        correspDesc._localEntityType, correspDesc._localGeoType, correspDesc._remoteEntityType,
        correspDesc._remoteGeoType, &toReturn[0] );

    return toReturn;
};

#endif
