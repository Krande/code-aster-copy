/**
 * @file ContactPairing.h
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

#pragma once

#include "Contact/ContactNew.h"
#include "Contact/ContactParameters.h"
#include "Contact/ContactZone.h"
#include "DataFields/FieldOnNodes.h"
#include "DataStructures/DataStructure.h"

class ContactPairing : public DataStructure {

  private:
    /** @brief new coordinates */
    MeshCoordinatesFieldPtr _newCoordinates;
    /** @brief Contact definition  */
    ContactNewPtr _contDefi;
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
    /** @brief Map between pair and zone */
    std::map< ASTERINTEGER, ASTERINTEGER > _pair2Zone;
    // FED after pairing
    FiniteElementDescriptorPtr _fed;

  public:
    typedef std::vector< std::pair< ASTERINTEGER, ASTERINTEGER > > VectorLongPairs;

    /** @brief Mesh constructor */
    ContactPairing( const std::string name, const ContactNewPtr cont );

    /** @brief constructor */
    ContactPairing( const ContactNewPtr cont )
        : ContactPairing( ResultNaming::getNewResultName(), cont ){};

    /** @brief Get coordinates */
    MeshCoordinatesFieldPtr getCoordinates() const { return _newCoordinates; }

    BaseMeshPtr getMesh() const { return _mesh; };

    /** @brief Update coordinates */
    void updateCoordinates( FieldOnNodesRealPtr &disp ) {
        *_newCoordinates = *( _mesh->getCoordinates() ) + *disp;
    };

    /** @brief Update coordinates */
    void setCoordinates( MeshCoordinatesFieldPtr &coor ) { _newCoordinates = coor; };

    /** @brief compute pairing quantities of zone i */
    ASTERBOOL computeZone( ASTERINTEGER i );

    ASTERBOOL compute();

    /** @brief clearZone all pairing quantities of zone i */
    void clearZone( ASTERINTEGER i );

    /** @brief clearZone all pairing quantities for all zones */
    void clear() {
        for ( auto i = 0; i < _contDefi->getNumberOfContactZones(); i++ ) {
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
    VectorLongPairs getListOfPairs() const;

    /** @brief get list of pairs of zone associated with zone zone_index
     *  @return vector of pairs of type std::pair
     */
    VectorLongPairs getListOfPairsOfZone( ASTERINTEGER zone_index ) const;

    /** @brief get slave intersection points of zone zone_index
     *   @return vector of slave intersection points of size 16*number of pairs
     **/
    std::vector< VectorReal > getSlaveIntersectionPoints( ASTERINTEGER zone_index ) const;

    /** @brief Get vector of number of intersection points  */
    std::vector< VectorLong > getNumberOfIntersectionPoints() const {
        return _nbIntersectionPoints;
    };

    /** @brief Get vector of slave intersection points */
    std::vector< VectorReal > getSlaveIntersectionPoints() const {
        return _slaveIntersectionPoints;
    };

    /** @brief Build Finite Element Descriptor from pairing */
    void buildFiniteElementDescriptor();

    /** @brief Get Finite Element Descriptor from pairing */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _fed; };

    /** @brief Get map */
    std::map< ASTERINTEGER, ASTERINTEGER > pairsToZones() const { return _pair2Zone; };
};

typedef std::shared_ptr< ContactPairing > ContactPairingPtr;