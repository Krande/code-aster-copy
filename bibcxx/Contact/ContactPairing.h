#ifndef PAIRING_H_
#define PAIRING_H_

/**
 * @file ContactPairing.h
 * @brief Fichier entete de la class Contact
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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

#include "Contact/ContactParameters.h"
#include "Contact/ContactZone.h"
#include "DataFields/FieldOnNodes.h"
#include "DataStructures/DataStructure.h"

class ContactPairing : public DataStructure {

  private:
    /** @brief new coordinates */
    MeshCoordinatesFieldPtr _newCoordinates;
    /** @brief vector of zones  */
    std::vector< ContactZonePtr > _zones;
    /** @brief Mesh */
    BaseMeshPtr _mesh;
    /** @brief Vector of number of pairs */
    VectorLong _nbPairs;
    /** @brief Vector of pairs */
    std::vector< VectorLong > _listOfPairs;
    /** @brief Vector of number of intersection points  */
    std::vector< VectorLong > _nbIntersectionPoints;
    /** @brief Vector of slave intersection points */
    std::vector< VectorReal > _slaveIntersectionPoints;
    /** @brief Vector of master intersection points */
    std::vector< VectorReal > _masterIntersectionPoints;
    /** @brief Gauss points */
    std::vector< VectorReal > _quadraturePoints;

  public:
    typedef std::vector< std::pair< ASTERINTEGER, ASTERINTEGER > > VectorLongPairs;

    /** @brief Mesh constructor */
    ContactPairing( const std::string name, const std::vector< ContactZonePtr > zones,
                    const BaseMeshPtr mesh );

    /** @brief Mesh constructor */
    ContactPairing( const std::vector< ContactZonePtr > zones, const BaseMeshPtr mesh )
        : ContactPairing( ResultNaming::getNewResultName(), zones, mesh ){};

    /** @brief Mesh getter */
    MeshCoordinatesFieldPtr getCoordinates() const { return _newCoordinates; }

    /** @brief zones  getter */
    std::vector< ContactZonePtr > getContactZones() const { return _zones; }

    /** @brief Update coordinates */
    void updateCoordinates( FieldOnNodesRealPtr &disp ) {
        *_newCoordinates = *( _mesh->getCoordinates() ) + *disp;
    };

    /** @brief compute pairing quantities of zone i */
    ASTERBOOL computeZone( ASTERINTEGER i );

    ASTERBOOL compute();

    /** @brief get zone of index zone_index */
    ContactZonePtr getContactZone( ASTERINTEGER zone_index ) { return _zones[zone_index]; }

    /** @brief clearZone all pairing quantities of zone i */
    void clearZone( ASTERINTEGER i );

    /** @brief clearZone all pairing quantities for all zones */
    void clear() {
        for ( auto i = 0; i < _zones.size(); i++ ) {
            clearZone( i );
        }
    };

    /** @brief get number of all pairs  */
    ASTERINTEGER getNumberOfPairs() const {

        return std::accumulate( _nbPairs.begin(), _nbPairs.end(), 0 );
    };

    /** @brief get number of pairs of zone zone_index  */
    ASTERINTEGER getNumberOfPairsOfZone( ASTERINTEGER zone_index ) const {
        return _nbPairs[zone_index];
    }

    /** @brief get list of pairs of zone associated with zone zone_index
     *  @return vector of pairs of type std::pair
     */
    VectorLongPairs getListOfPairsOfZone( ASTERINTEGER zone_index ) const;

    /** @brief get slave intersection points of zone zone_index
     *   @return vector of slave intersection points of size 16*number of pairs
     **/
    VectorReal getSlaveIntersectionPoints( ASTERINTEGER zone_index ) const {
        return _slaveIntersectionPoints[zone_index];
    }

    /** @brief get master intersection points of zone zone_index
     *  @return vector of master intersection points of size 16*number of pairs
     **/
    VectorReal getMasterIntersectionPoints( ASTERINTEGER zone_index ) const {
        return _masterIntersectionPoints[zone_index];
    }

    /** @brief get gauss points of zone zone_index
     * @return vector Gauss quadrature points of size 72*number of pairs
     **/
    VectorReal getQuadraturePoints( ASTERINTEGER zone_index ) const {
        return _quadraturePoints[zone_index];
    }
};

typedef std::shared_ptr< ContactPairing > ContactPairingPtr;

#endif /* PAIRING_H_ */
