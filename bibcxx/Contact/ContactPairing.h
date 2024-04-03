/**
 * @file ContactPairing.h
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

#pragma once

#include "Contact/ContactNew.h"
#include "Contact/ContactParameter.h"
#include "Contact/ContactZone.h"
#include "DataFields/FieldOnNodes.h"
#include "DataStructures/DataStructure.h"

// Description of a virtual contact cell
struct contCellType {
    std::string contCellType; // Geometric support
    int nbNode;               // Number of nodes of geometric support
    std::string slavCellType; // Geometric support for slave side
    std::string mastCellType; // Geometric support for master side
    std::string contElemType; // Finite element for contact
    std::string fricElemType; // Finite element for firrction
};

// For ALGO_CONT="LAGRANGIEN"
const unsigned int contLagrType = 33;
const struct contCellType contCellLagr[contLagrType] = {
    { "SEG22", 4, "SEG2", "SEG2", "CMS2S2", "FMS2S2" },
    { "SEG33", 6, "SEG3", "SEG3", "CMS3S3", "FMS3S3" },
    { "SEG23", 5, "SEG2", "SEG3", "CMS2S3", "FMS2S3" },
    { "SEG32", 5, "SEG3", "SEG2", "CMS3S2", "FMS3S2" },
    { "TRIA33", 6, "TRIA3", "TRIA3", "CMT3T3", "FMT3T3" },
    { "TR3TR6", 9, "TRIA3", "TRIA6", "CMT3T6", "FMT3T6" },
    { "TR6TR3", 9, "TRIA6", "TRIA3", "CMT6T3", "FMT6T3" },
    { "TRIA66", 12, "TRIA6", "TRIA6", "CMT6T6", "FMT6T6" },
    { "QUAD44", 8, "QUAD4", "QUAD4", "CMQ4Q4", "FMQ4Q4" },
    { "QU4QU8", 12, "QUAD4", "QUAD8", "CMQ4Q8", "FMQ4Q8" },
    { "QU8QU4", 12, "QUAD8", "QUAD4", "CMQ8Q4", "FMQ8Q4" },
    { "QUAD88", 16, "QUAD8", "QUAD8", "CMQ8Q8", "FMQ8Q8" },
    { "QU4TR3", 7, "QUAD4", "TRIA3", "CMQ4T3", "FMQ4T3" },
    { "TR3QU4", 7, "TRIA3", "QUAD4", "CMT3Q4", "FMT3Q4" },
    { "TR6QU4", 10, "TRIA6", "QUAD4", "CMT6Q4", "FMT6Q4" },
    { "QU4TR6", 10, "QUAD4", "TRIA6", "CMQ4T6", "FMQ4T6" },
    { "TR6QU8", 14, "TRIA6", "QUAD8", "CMT6Q8", "FMT6Q8" },
    { "QU8TR6", 14, "QUAD8", "TRIA6", "CMQ8T6", "FMQ8T6" },
    { "TR6QU9", 15, "TRIA6", "QUAD9", "CMT6Q9", "FMT6Q9" },
    { "QU9TR6", 15, "QUAD9", "TRIA6", "CMQ9T6", "FMQ9T6" },
    { "QU8TR3", 11, "QUAD8", "TRIA3", "CMQ8T3", "FMQ8T3" },
    { "TR3QU8", 11, "TRIA3", "QUAD8", "CMT3Q8", "FMT3Q8" },
    { "QU8QU9", 17, "QUAD8", "QUAD9", "CMQ8Q9", "FMQ8Q9" },
    { "QU9QU8", 17, "QUAD9", "QUAD8", "CMQ9Q8", "FMQ9Q8" },
    { "QU9QU4", 13, "QUAD9", "QUAD4", "CMQ9Q4", "FMQ9Q4" },
    { "QU4QU9", 13, "QUAD4", "QUAD9", "CMQ4Q9", "FMQ4Q9" },
    { "QU9TR3", 12, "QUAD9", "TRIA3", "CMQ9T3", "FMQ9T3" },
    { "TR3QU9", 12, "TRIA3", "QUAD9", "CMT3Q9", "FMT3Q9" },
    { "QUAD99", 18, "QUAD9", "QUAD9", "CMQ9Q9", "FMQ9Q9" },
    { "POI1", 1, "POI1", "LAG2", "CMP1L2", "FMP1L2" },
    { "POI1", 1, "POI1", "NOLAG2", "CMP1N2", "FMP1N2" },
    { "POI1", 1, "POI1", "LAG3", "CMP1L3", "FMP1L3" },
    { "POI1", 1, "POI1", "NOLAG3", "CMP1N3", "FMP1N3" },
};

// For ALGO_CONT="NITSCHE"
const unsigned int contNitsType = 10;
const struct contCellType contCellNits[contNitsType] = {
    { "TR3SE2", 5, "TRIA3", "SEG2", "CNT3S2", "FNT3S2" },
    { "TR3SE3", 6, "TRIA3", "SEG3", "CNT3S3", "FNT3S3" },
    { "TR6SE2", 8, "TRIA6", "SEG2", "CNT6S2", "FNT6S2" },
    { "TR6SE3", 9, "TRIA6", "SEG3", "CNT6S3", "FNT6S3" },
    { "QU4SE2", 6, "QUAD4", "SEG2", "CNQ4S2", "FNQ4S2" },
    { "QU4SE3", 7, "QUAD4", "SEG3", "CNQ4S3", "FNQ4S3" },
    { "QU8SE2", 10, "QUAD8", "SEG2", "CNQ8S2", "FNQ8S2" },
    { "QU8SE3", 11, "QUAD8", "SEG3", "CNQ8S3", "FNQ8S3" },
    { "QU9SE2", 11, "QUAD9", "SEG2", "CNQ9S2", "FNQ9S2" },
    { "QU9SE3", 12, "QUAD9", "SEG3", "CNM9S3", "FNM9S3" },
};

class ContactPairing : public DataStructure {
    /** Datastructure for pairing */
  private:
    /** @brief Current coordinates of nodes */
    MeshCoordinatesFieldPtr _currentCoordinates;

    /** @brief Contact definition */
    ContactNewPtr _contDefi;

    /** @brief Mesh */
    BaseMeshPtr _mesh;

    /** @brief Vector of number of pairs */
    VectorLong _nbPairs;

    /** @brief Vector of pairs */
    std::vector< VectorLong > _listOfPairs;

    /** @brief Vector of number of intersection points  */
    std::vector< VectorLong > _nbIntersectionPoints;

    /** @brief Vector of coordinates for intersection points */
    std::vector< VectorReal > _slaveIntersectionPoints;

    /** @brief Map between pair and zone */
    MapLong _pair2Zone;

    /** @brief Finite element descriptor for virtual elements of contact */
    FiniteElementDescriptorPtr _fed;

  private:
    /** @brief Resize pairing quantities */
    void resizePairing( const int nbZoneCont );

    /** @brief Create virtual elements for contact */
    void createVirtualElemForContact( const ASTERLOGICAL lAxis, const int nbZoneCont,
                                      MapLong &contactElemType,
                                      const JeveuxCollectionLong meshConnectivity,
                                      std::vector< VectorLong > &listContElem,
                                      std::vector< VectorPairLong > &listContType,
                                      ASTERINTEGER &iContPair, SetLong &slaveNodePaired,
                                      SetLong &slaveCellPaired );

    /** @brief Create virtual elements for orphelan nodes */
    void createVirtualElemForOrphelanNodes( const ASTERLOGICAL lAxis, const int nbZoneCont,
                                            MapLong &contactElemType,
                                            const JeveuxCollectionLong meshConnectivity,
                                            std::vector< VectorLong > &listContElem,
                                            std::vector< VectorPairLong > &listContType,
                                            ASTERINTEGER &iContPair, SetLong &slaveNodePaired,
                                            SetLong &slaveCellPaired );

    /** @brief Get index of contact cell */
    ASTERINTEGER getContCellIndx( const ContactAlgo contAlgo, std::string slavCellTypeName,
                                  std::string mastCellTypeName );

    /** @brief Get type of contact cell */
    ASTERINTEGER getContCellType( const ContactAlgo contAlgo, const ASTERINTEGER cellIndx,
                                  const bool lAxis, const bool lFric );

  public:
    using VectorLongPairs = std::vector< std::pair< ASTERINTEGER, ASTERINTEGER > >;

    /** @brief Constructor with given name */
    ContactPairing( const std::string name, const ContactNewPtr cont );

    /** @brief Constructor with automatic name */
    ContactPairing( const ContactNewPtr cont )
        : ContactPairing( ResultNaming::getNewResultName(), cont ) {};

    /** @brief Get coordinates */
    MeshCoordinatesFieldPtr getCoordinates() const { return _currentCoordinates; }

    /** @brief Get mesh */
    BaseMeshPtr getMesh() const { return _mesh; };

    /** @brief Update coordinates */
    void updateCoordinates( FieldOnNodesRealPtr &disp ) {
        *_currentCoordinates = *( _mesh->getCoordinates() ) + *disp;
    };

    /** @brief Set coordinates */
    void setCoordinates( MeshCoordinatesFieldPtr &coor ) { _currentCoordinates = coor; };

    /** @brief Compute pairing quantities of zone */
    ASTERBOOL computeZone( ASTERINTEGER indexZone );

    /** @brief Compute pairing quantities of all zones */
    ASTERBOOL compute();

    /** @brief Clear pairing quantities of zone */
    void clearZone( ASTERINTEGER indexZone );

    /** @brief Clear pairing quantities for all zones */
    void clear() {
        for ( auto indexZone = 0; indexZone < _contDefi->getNumberOfContactZones(); indexZone++ ) {
            clearZone( indexZone );
        }
    };

    /** @brief Get number of all pairs  */
    ASTERINTEGER getNumberOfPairs() const {
        return std::accumulate( _nbPairs.begin(), _nbPairs.end(), 0 );
    };

    /** @brief Get number of pairs of zone  */
    ASTERINTEGER getNumberOfPairsOfZone( ASTERINTEGER indexZone ) const {
        return _nbPairs[indexZone];
    }

    /** @brief Get all list of pairs */
    VectorLongPairs getListOfPairs() const;

    /** @brief Get list of pairs of zone  */
    VectorLongPairs getListOfPairsOfZone( ASTERINTEGER indexZone ) const;

    /** @brief Get slave intersection points of zone */
    std::vector< VectorReal > getSlaveIntersectionPoints( ASTERINTEGER indexZone ) const;

    /** @brief Get vector of number of intersection points on all zones */
    std::vector< VectorLong > getNumberOfIntersectionPoints() const {
        return _nbIntersectionPoints;
    };

    /** @brief Get vector of slave intersection points on all zones */
    std::vector< VectorReal > getSlaveIntersectionPoints() const {
        return _slaveIntersectionPoints;
    };

    /** @brief Build Finite Element Descriptor from pairing */
    void buildFiniteElementDescriptor();

    /** @brief Get Finite Element Descriptor from pairing */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _fed; };

    /** @brief Get map between pair and zone */
    MapLong pairsToZones() const { return _pair2Zone; };
};

using ContactPairingPtr = std::shared_ptr< ContactPairing >;
