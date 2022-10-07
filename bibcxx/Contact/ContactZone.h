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
    /** @brief  Check direction of normal */
    bool _checkNormal;
    /** @brief  Smoothing of normal */
    bool _smoothing;
    /** @brief  List of master Cells */
    VectorLong _masterCells;
    /** @brief List of slave cells */
    VectorLong _slaveCells;
    /** @brief  List of master nodes */
    VectorLong _masterNodes;
    /** @brief List of slave nodes */
    VectorLong _slaveNodes;
    /** @brief excluded elements of slave side for LAGRANGIEN */
    VectorLong _slaveCellsExcluded;
    /** @brief  Master inverse connectivity */
    JeveuxCollectionLong _masterInverseConnectivity;
    /** @brief  Slave inverse connectivity */
    JeveuxCollectionLong _slaveInverseConnectivity;
    /** @brief  Master cells neighbors */
    JeveuxCollectionLong _masterNeighbors;
    /** @brief  slave cells neighbors */
    JeveuxCollectionLong _slaveNeighbors;
    /** @brief Map between slave surfaciv and volumic cell */
    std::map< ASTERINTEGER, ASTERINTEGER > _slavSurf2Volu;
    /** @brief name of slave side */
    std::string _slaveGrp;
    /** @brief name of master side */
    std::string _masterGrp;

    /**
     * @brief Construct the inverse connectivity
     */
    ASTERBOOL buildInverseConnectivity();
    /**
     * @brief construct master/slave cells neighbors
     */
    ASTERBOOL buildCellsNeighbors();

    /**
     * @brief construct surface->volume slave cell
     */
    void buildSlaveCellsVolu();

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

    void setSlaveGroupOfCells( const std::string &slave );

    void setMasterGroupOfCells( const std::string &master );

    void setExcludedSlaveGroupOfCells( const VectorString &excluded_slave );

    void setExcludedSlaveGroupOfNodes( const VectorString &excluded_slave );

    VectorLong getExcludedSlaveCells() const { return _slaveCellsExcluded; };

    void checkNormals( const bool &checkNormal ) { _checkNormal = checkNormal; }

    bool checkNormals() const { return _checkNormal; }

    /**
     * @brief get master nodes
     */
    const VectorLong &getMasterNodes() const { return _masterNodes; };

    VectorLong &getMasterNodes() {
        return const_cast< VectorLong & >( std::as_const( *this ).getMasterNodes() );
    }

    VectorLong getSlaveNodes() const { return _slaveNodes; }

    /**
     * @brief get master cells
     */
    const VectorLong &getMasterCells() const { return _masterCells; };

    VectorLong &getMasterCells() {
        return const_cast< VectorLong & >( std::as_const( *this ).getMasterCells() );
    }

    /**
     * @brief get slave cells
     */
    VectorLong getSlaveCells() const { return _slaveCells; }

    VectorLong getMasterCellsFromNode( const ASTERINTEGER &i ) const;

    VectorLong getSlaveCellsFromNode( const ASTERINTEGER &i ) const;

    VectorLong getMasterCellNeighbors( const ASTERINTEGER &i ) const;

    VectorLong getSlaveCellNeighbors( const ASTERINTEGER &i ) const;

    auto getSlaveCellsSurfToVolu() const { return _slavSurf2Volu; };

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

    void enableFriction( const bool &friction ) { _fricParam->enableFriction( friction ); };

    bool hasFriction() const { return _fricParam->hasFriction(); };

    void enableSmoothing( const bool &smoothing ) { _smoothing = smoothing; };

    bool hasSmoothing() const { return _smoothing; };
};

/**
 * @typedef ContactZonePtr
 * @brief Pointeur intelligent vers un ContactZone
 */
typedef std::shared_ptr< ContactZone > ContactZonePtr;

#endif /* CONTACT_ZONE_H_ */
