#ifndef CONTACT_ZONE_H_
#define CONTACT_ZONE_H_

/**
 * @file ContactZone.h
 * @brief Fichier entete de la class ContactZone
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

#include "astercxx.h"

#include "Contact/ContactEnum.h"
#include "Contact/ContactParameters.h"
#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"

class ContactZone : public DataStructure {
  private:
    /** @brief Modele */
    ModelPtr _model;
    /** @brief Level of verbosity */
    ASTERINTEGER _verbosity;
    /** @brief Parameter for contact only */
    ContactParameterPtr _contParam;
    /** @brief Parameter for friction only */
    FrictionParameterPtr _fricParam;
    /** @brief Parameter for pairing only */
    PairingParameterPtr _pairParam;
    /** @brief Slave side */
    std::string _slave;
    /** @brief Master side */
    std::string _master;
    /** @brief excluded elements of slave side for LAGRANGIEN */
    VectorString _excluded_slave;
    /** @brief  Check direction of normal */
    bool _checkNormal;
    /** @brief  List of master Cells */
    VectorLong _masterCells;
    /** @brief List of slave cells */
    VectorLong _slaveCells;
    /** @brief  Master inverse connectivity */
    JeveuxCollectionLong _masterInverseConnectivity;
    /** @brief  Slave inverse connectivity */
    JeveuxCollectionLong _slaveInverseConnectivity;
    /** @brief  Master Nodes */
    VectorLong _masterNodes;
    /** @brief  Master cells neighbors */
    JeveuxCollectionLong _masterNeighbors;
    /** @brief  slave cells neighbors */
    JeveuxCollectionLong _slaveNeighbors;

  public:
    /**
     * @typedef ContactZonePt
     * @brief Pointeur intelligent vers un ContactZone
     */
    typedef std::shared_ptr< ContactZone > ContactZonePtr;
    /**
     * @brief Constructeur
     */
    ContactZone() = delete;

    /**
     * @brief Constructeur
     */
    ContactZone( const std::string name, const ModelPtr model );

    /**
     * @brief Constructeur
     */
    ContactZone( const ModelPtr model ) : ContactZone( ResultNaming::getNewResultName(), model ){};

    ModelPtr getModel() const { return _model; }

    BaseMeshPtr getMesh() const { return _model->getMesh(); }

    void setVerbosity( const ASTERINTEGER &level ) { _verbosity = level; }

    ASTERINTEGER getVerbosity() const { return _verbosity; }

    bool build();

    ContactParameterPtr getContactParameter() const { return _contParam; };

    FrictionParameterPtr getFrictionParameter() const { return _fricParam; };

    PairingParameterPtr getPairingParameter() const { return _pairParam; };

    void setContactParameter( const ContactParameterPtr contParam ) { _contParam = contParam; };

    void setFrictionParameter( const FrictionParameterPtr fricParam ) { _fricParam = fricParam; };

    void setPairingParameter( const PairingParameterPtr pairParam ) { _pairParam = pairParam; };

    void setSlaveGroupOfCells( const std::string &slave ) {
        if ( getMesh()->hasGroupOfCells( slave ) ) {
            _slave = slave;
        } else {
            throw std::runtime_error( "The given group " + slave + " doesn't exist in mesh" );
        }
    };

    std::string getSlaveGroupOfCells() const { return _slave; };

    void setMasterGroupOfCells( const std::string &master ) {
        if ( getMesh()->hasGroupOfCells( master ) ) {
            _master = master;
        } else {
            throw std::runtime_error( "The given group " + master + " doesn't exist in mesh" );
        }
    };

    std::string getMasterGroupOfCells() const { return _master; };

    void setExcludedSlaveGroupOfCells( const VectorString &excluded_slave ) {
        for ( long i = 0; i < excluded_slave.size(); i++ ) {
            if ( !( getMesh()->hasGroupOfCells( excluded_slave[i] ) ) ) {
                throw std::runtime_error( "The group " + excluded_slave[0] +
                                          " doesn't exist in mesh" );
            }
        }
        _excluded_slave = excluded_slave;
    };

    VectorString getExcludedSlaveGroupOfCells() const { return _excluded_slave; };

    void checkNormals( const bool &checkNormal ) { _checkNormal = checkNormal; }

    bool checkNormals() const { return _checkNormal; }

    /**
     * @brief get master nodes
     */
    const VectorLong &getMasterNodes() const { return _masterCells; };

    VectorLong &getMasterNodes() {
        return const_cast< VectorLong & >( std::as_const( *this ).getMasterNodes() );
    }

    /**
     * @brief get master cells
     */
    const VectorLong &getMasterCells() const { return _masterCells; };

    VectorLong &getMasterCells() {
        return const_cast< VectorLong & >( std::as_const( *this ).getMasterCells() );
    }

    /**
     * @brief  update list of master cells
     */
    void updateMasterCells() {
        if ( _masterCells.empty() )
            _masterCells = getMesh()->getCells( _master );
    }

    /**
     * @brief get slave cells
     */
    const VectorLong &getSlaveCells() const { return _slaveCells; };

    VectorLong &getSlaveCells() {
        return const_cast< VectorLong & >( std::as_const( *this ).getSlaveCells() );
    }

    /**
     * @brief update list of slave cells
     */
    void updateSlaveCells() {
        if ( _slaveCells.empty() )
            _slaveCells = getMesh()->getCells( _slave );
    }

    VectorLong getMasterCellsFromNode( const int &i ) const {
        auto vct = _masterInverseConnectivity->getObject( i ).toVector();
        std::transform(
            vct.begin(), vct.end(), vct.begin(),
            [this]( ASTERINTEGER k ) -> ASTERINTEGER { return k > 0 ? _masterCells[k - 1] : 0; } );
        return vct;
    }

    VectorLong getSlaveCellsFromNode( const int &i ) const {
        auto vct = _slaveInverseConnectivity->getObject( i ).toVector();
        std::transform(
            vct.begin(), vct.end(), vct.begin(),
            [this]( ASTERINTEGER k ) -> ASTERINTEGER { return k > 0 ? _slaveCells[k - 1] : 0; } );
        return vct;
    }

    VectorLong getMasterCellNeighbors( const int &i ) const {
        ASTERINTEGER ind_min = *std::min_element( _masterCells.begin(), _masterCells.end() );
        ASTERINTEGER ind_max = *std::max_element( _masterCells.begin(), _masterCells.end() );

        if ( i < ind_min || i > ind_max )
            throw std::out_of_range( " the master cell's number should be"
                                     " between " +
                                     std::to_string( ind_min ) + " and " +
                                     std::to_string( ind_max ) );

        auto vct = _masterNeighbors->getObject( i - ind_min + 1 ).toVector();
        vct.erase(
            std::remove_if( vct.begin(), vct.end(), []( ASTERINTEGER &i ) { return i == 0; } ),
            vct.end() );
        return vct;
    }

    VectorLong getSlaveCellNeighbors( const int &i ) const {

        ASTERINTEGER ind_min = *std::min_element( _slaveCells.begin(), _slaveCells.end() );
        ASTERINTEGER ind_max = *std::max_element( _slaveCells.begin(), _slaveCells.end() );

        if ( i < ind_min || i > ind_max )
            throw std::out_of_range( " the slave cell's number should be"
                                     " between " +
                                     std::to_string( ind_min ) + " and " +
                                     std::to_string( ind_max ) );

        auto vct = _slaveNeighbors->getObject( i - ind_min + 1 ).toVector();
        vct.erase(
            std::remove_if( vct.begin(), vct.end(), []( ASTERINTEGER &i ) { return i == 0; } ),
            vct.end() );
        return vct;
    }

    /**
     * @brief get master inverse connectivity as JeVeuxCollection
     */
    JeveuxCollectionLong getMasterInverseConnectivity() const { return _masterInverseConnectivity; }
    /**
     * @brief get slave inverse connectivity as JeVeuxCollection
     */
    JeveuxCollectionLong getSlaveInverseConnectivity() const { return _slaveInverseConnectivity; }

    /**
     * @brief get master neighbors
     */
    JeveuxCollectionLong getMasterNeighbors() const { return _masterNeighbors; }
    /**
     * @brief get slave neighbors
     */
    JeveuxCollectionLong getSlaveNeighbors() const { return _slaveNeighbors; }
    /**
     * @brief Construct the inverse connectivity
     */
    ASTERBOOL buildInverseConnectivity();
    /**
     * @brief construct master/slave cells neighbors
     */
    ASTERBOOL buildCellsNeighbors();
};

/**
 * @typedef ContactZonePtr
 * @brief Pointeur intelligent vers un ContactZone
 */
typedef std::shared_ptr< ContactZone > ContactZonePtr;

#endif /* CONTACT_ZONE_H_ */
