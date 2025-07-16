#ifndef MESH_PAIRING_H_
#define MESH_PAIRING_H_

/**
 * @file MeshPairing.h
 * @brief Header of MeshPairing class
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

#include "astercxx.h"

#include "DataFields/FieldOnNodes.h"
#include "DataFields/MeshCoordinatesField.h"
#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/Mesh.h"
#include "Meshes/MeshEnum.h"
#include "Supervis/ResultNaming.h"

class MeshPairing : public DSWithCppPickling {
  private:
    /** @brief Mesh */
    BaseMeshPtr _mesh;
    /** @brief Current coordinates of nodes */
    MeshCoordinatesFieldPtr _currentCoordinates;

    /** @brief Flag when zone is defined */
    bool _zoneHaveBeenDefined;
    /** @brief Name of group for slave side */
    std::string _slaveCellsGroup;
    /** @brief Name of group for master side */
    std::string _masterCellsGroup;
    /** @brief List of master cells */
    VectorLong _masterCells;
    /** @brief List of slave cells */
    VectorLong _slaveCells;
    /** @brief List of master nodes */
    VectorLong _masterNodes;
    /** @brief List of slave nodes */
    VectorLong _slaveNodes;
    /** @brief List of excluded slave cells */
    VectorLong _slaveCellsExcluded;
    /** @brief List of groups of excluded slave cells */
    VectorString _excludedSlaveCells;
    /** @brief List of groups of excluded slave nodes */
    VectorString _excludedSlaveNodes;

    /** @brief  Master inverse connectivity */
    JeveuxCollectionLong _masterInverseConnectivity;
    /** @brief  Slave inverse connectivity */
    JeveuxCollectionLong _slaveInverseConnectivity;
    /** @brief  Master cells neighbors */
    JeveuxCollectionLong _masterNeighbors;
    /** @brief  Slave cells neighbors */
    JeveuxCollectionLong _slaveNeighbors;

    /** @brief  Output JEVEUX objects from FORTRAN pairing */
    JeveuxVectorLong _jvPairs;
    JeveuxVectorLong _jvNbInterPoints;
    JeveuxVectorReal _jvInterSlavePoints;

    /** @brief Map between slave surfacic and volumic cell */
    MapLong _slavSurf2Volu;

    /** @brief Method of pairing */
    PairingMethod _method;

    /** @brief Level of verbosity */
    ASTERINTEGER _verbosity;

    /** @brief Number of pairs */
    ASTERINTEGER _nbPairs;

    /** @brief Vector of pairs */
    VectorLong _pairs;

    /** @brief Vector of number of intersection points  */
    VectorLong _nbPoinInte;

    /** @brief Vector of coordinates for intersection points in slave parameteric space */
    VectorReal _poinInteSlav;

  private:
    /** @brief Get name of JEVEUX object for master inverse connectivity */
    std::string getMasterInverseConnectivityName() const {
        return ljust( getName() + ".CM", 24, ' ' );
    }

    /** @brief Get name of JEVEUX object for slave inverse connectivity */
    std::string getSlaveInverseConnectivityName() const {
        return ljust( getName() + ".CS", 24, ' ' );
    }

    /** @brief Get name of JEVEUX object for master neighbors */
    std::string getMasterNeighName() const { return ljust( getName() + ".MN", 24, ' ' ); }

    /** @brief Get name of JEVEUX object for slave neighbors */
    std::string getSlaveNeighName() const { return ljust( getName() + ".SN", 24, ' ' ); }

    /** @brief Get name of JEVEUX object for outputs*/
    std::string getBasename() const { return ljust( getName(), 8, ' ' ); }

    /** @brief Get name of JEVEUX object for pairs */
    std::string getPairsName() const { return ljust( getBasename() + ".LISTPAIRS", 24, ' ' ); }

    /** @brief Get name of JEVEUX object for number of intersection points */
    std::string getNbInterName() const { return ljust( getBasename() + ".NBPOIN", 24, ' ' ); }

    /** @brief Get name of JEVEUX object for coordinates of intersection points */
    std::string getCoorInterName() const { return ljust( getBasename() + ".INTERSLPTS", 24, ' ' ); }

    /** @brief Construct the inverse connectivity */
    ASTERBOOL buildInverseConnectivity();

    /** @brief Construct master/slave cells neighbors */
    ASTERBOOL buildCellsNeighbors();

    /** @brief Construct surface to volume slave cell mapping */
    ASTERBOOL buildSlaveCellsVolu();

    /** @brief Surface defined */
    ASTERBOOL surfacesHasBeenDefined();

  public:
    using MeshPairingPtr = std::shared_ptr< MeshPairing >;

    /** @brief Constructor with given name */
    MeshPairing( const std::string name );

    /** @brief Constructor with automatic name */
    MeshPairing() : MeshPairing( ResultNaming::getNewResultName() ) {};

    /** @brief restricted constructor (Set) and method (Get) to support pickling */
    MeshPairing( const py::tuple &tup );
    py::tuple _getState() const;

    /** @brief Initializations of datastructures defining pairing */
    bool build();

    /** @brief Get mesh */
    BaseMeshPtr getMesh() const { return _mesh; };
    void setMesh( const BaseMeshPtr & );

    /** @brief Set verbosity */
    void setVerbosity( const ASTERINTEGER &level ) { _verbosity = level; }

    /** @brief Get verbosity */
    ASTERINTEGER getVerbosity() const { return _verbosity; }

    /** @brief Get coordinates */
    MeshCoordinatesFieldPtr getCoordinates() const { return _currentCoordinates; }

    // /** @brief Set coordinates */
    void setCoordinates( const MeshCoordinatesFieldPtr &coor ) { _currentCoordinates = coor; };

    /** @brief Set pair */
    void setPair( const std::string &groupNameSlav, const std::string &groupNameMast );

    /** @brief Set slave cells */
    void setSlaveGroupOfCells( const std::string &groupName );

    /** @brief Set master cells */
    void setMasterGroupOfCells( const std::string &groupName );

    /** @brief Set excluded groups of slave cells */
    void setExcludedSlaveGroupOfCells( const VectorString &groupsName );

    /** @brief Set excluded groups of slave nodes */
    void setExcludedSlaveGroupOfNodes( const VectorString &groupsName );

    /** @brief Get name of group for slave cells */
    const std::string &getSlaveGroupOfCells() const { return _slaveCellsGroup; }

    /** @brief Get name of group for master cells */
    const std::string &getMasterGroupOfCells() const { return _masterCellsGroup; }

    /** @brief Get master cells*/
    const VectorLong &getMasterCells() const { return _masterCells; };

    /** @brief Get slave cells */
    const VectorLong &getSlaveCells() const { return _slaveCells; };

    /** @brief Get master nodes*/
    const VectorLong &getMasterNodes() const { return _masterNodes; };
    VectorLong &getMasterNodes() {
        return const_cast< VectorLong & >( std::as_const( *this ).getMasterNodes() );
    };

    /** @brief Get slave nodes*/
    const VectorLong &getSlaveNodes() const { return _slaveNodes; };
    VectorLong &getSlaveNodes() {
        return const_cast< VectorLong & >( std::as_const( *this ).getSlaveNodes() );
    };

    /** @brief Get master cells from a node (inverse connectivity) */
    VectorLong getMasterCellsFromNode( const ASTERINTEGER &nodeIndex ) const;

    /** @brief Get slave cells from a node (inverse connectivity) */
    VectorLong getSlaveCellsFromNode( const ASTERINTEGER &nodeIndex ) const;

    /** @brief Get neighbors of master cell */
    VectorLong getMasterCellNeighbors( const ASTERINTEGER &cellIndex ) const;

    /** @brief Get neighbors of slave cell */
    VectorLong getSlaveCellNeighbors( const ASTERINTEGER &cellIndex ) const;

    /** @brief Get volume slave cells linked to all surfacic slave cells */
    MapLong getSlaveCellsSurfToVolu() const { return _slavSurf2Volu; };

    /** @brief Get volume slave cells linked to one surfacic slave cells */
    ASTERINTEGER getSlaveCellSurfToVolu( const ASTERINTEGER &cellIndex ) const {
        return _slavSurf2Volu.at( cellIndex );
    };

    /** @brief Main subroutine for pairing */
    ASTERBOOL compute( ASTERDOUBLE &dist_pairing, ASTERDOUBLE &pair_tole );

    /** @brief Clear pairing result */
    void clearResult();

    /** @brief Get number of pairs  */
    ASTERBOOL hasPairs() const { return ( _nbPairs > 0 ); };

    /** @brief Get number of pairs  */
    ASTERINTEGER getNumberOfPairs() const { return _nbPairs; };

    /** @brief Get all list of pairs */
    VectorPairLong getListOfPairs() const;

    /** @brief Get pair */
    MapLong getPair( const ASTERINTEGER &iPair );

    /** @brief Get number of intersection points on all pairs */
    VectorLong getNumberOfIntersectionPoints() const;

    /** @brief Get number of intersection points of given pair  */
    ASTERINTEGER getNumberOfIntersectionPoints( const ASTERINTEGER &indexPair ) const;

    /** @brief Get intersection points of given pair  */
    VectorOfVectorsReal
    getIntersectionPoints( const ASTERINTEGER &indexPair,
                           const CoordinatesSpace = CoordinatesSpace::Global ) const;

    /** @brief Get area of intersection of given pair  */
    ASTERDOUBLE getIntersectionArea( const ASTERINTEGER &indexPair ) const;

    /** @brief Get number of quadrature points of given pair  */
    ASTERINTEGER getNumberOfQuadraturePoints( const ASTERINTEGER &indexPair ) const;

    /** @brief Get quadrature points of given pair  */
    VectorOfVectorsReal getQuadraturePoints( const ASTERINTEGER &indexPair ) const;

    /** @brief Check for common nodes between slave and master side */
    ASTERBOOL hasCommonNodes() const;

    /** @brief Check orientation of normals */
    void checkNormals( const ModelPtr _model ) const;

    /** @brief Set method */
    void setMethod( const PairingMethod &method ) { _method = method; };
};

using MeshPairingPtr = std::shared_ptr< MeshPairing >;

#endif
