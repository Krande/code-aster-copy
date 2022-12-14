#ifndef BASEMESH_H_
#define BASEMESH_H_

/**
 * @file BaseMesh.h
 * @brief Fichier entete de la classe BaseMesh
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

#include "DataFields/ListOfTables.h"
#include "DataFields/MeshCoordinatesField.h"
#include "DataStructures/DataStructure.h"
#include "MemoryManager/NamesMap.h"
#include "Meshes/MeshExplorer.h"
#include "Utilities/GenericEnum.h"

/** @brief Forward declaration of ConstantFieldOnCells */
template < class ValueType >
class ConstantFieldOnCells;
typedef ConstantFieldOnCells< ASTERDOUBLE > ConstantFieldOnCellsReal;
typedef std::shared_ptr< ConstantFieldOnCellsReal > ConstantFieldOnCellsRealPtr;

/**
 * @class BaseMesh
 * @brief This object is the base class for all meshes variants
 */
class BaseMesh : public DataStructure, public ListOfTables {
  public:
    typedef MeshExplorer< CellsIteratorFromConnectivity, const JeveuxCollectionLong &,
                          const JeveuxVectorLong & >
        ConnectivityMeshExplorer;

  protected:
    typedef JeveuxCollection< ASTERINTEGER, NamesMapChar24 > JeveuxCollectionLongNamePtr;
    /** @brief Objet Jeveux '.DIME' */
    JeveuxVectorLong _dimensionInformations;
    /** @brief Pointeur de nom Jeveux '.NOMNOE' */
    NamesMapChar8 _nameOfNodes;
    /** @brief Champ aux noeuds '.COORDO' */
    MeshCoordinatesFieldPtr _coordinates;
    /** @brief Pointeur de nom Jeveux '.PTRNOMNOE' */
    NamesMapChar24 _nameOfGrpNodes;
    /** @brief Collection Jeveux '.GROUPENO' */
    JeveuxCollectionLongNamePtr _groupsOfNodes;
    /** @brief Collection Jeveux '.CONNEX' */
    JeveuxCollectionLong _connectivity;
    /** @brief Pointeur de nom Jeveux '.NOMMAIL' */
    NamesMapChar8 _nameOfCells;
    /** @brief Objet Jeveux '.TYPMAIL' */
    JeveuxVectorLong _cellsType;
    /** @brief Pointeur de nom Jeveux '.PTRNOMMAI' */
    NamesMapChar24 _nameOfGrpCells;
    /** @brief Objet Jeveux '.GROUPEMA' */
    JeveuxCollectionLongNamePtr _groupsOfCells;
    /** @brief jeveux vector '.ADAPTATION' */
    JeveuxVectorLong _adapt;
    /** @brief jeveux vector '.MAOR' */
    JeveuxVectorChar8 _oriMeshName;
    /** @brief jeveux vector '.CRMA' */
    JeveuxVectorLong _oriMeshCells;
    /** @brief jeveux vector '.CRNO' */
    JeveuxVectorLong _oriMeshNodes;
    /** @brief Collection jeveux '.PATCH' */
    JeveuxCollectionLong _patch;
    /** @brief jeveux vector '.CONOPA' */
    JeveuxVectorLong _nodePatchConnectivity;
    /** @brief jeveux vector '.COMAPA' */
    JeveuxVectorLong _cellPatchConnectivity;
    /** @brief jeveux vector '.PTRNOMPAT' */
    JeveuxVectorChar24 _namePatch;
    /** @brief jeveux vector '.NOMACR' */
    JeveuxVectorLong _superElementName;
    /** @brief jeveux vector '.PARA_R' */
    JeveuxVectorReal _superElementPara;
    /** @brief jeveux vector '.SUPMAIL' */
    JeveuxCollectionLong _superElements;
    /** @brief card '.ABSC_CURV' */
    ConstantFieldOnCellsRealPtr _curvAbsc;
    /** @brief Object to allow loop over connectivity */
    const ConnectivityMeshExplorer _explorer;

    /**
     * @brief Constructeur
     * @param name nom jeveux de l'objet
     * @param type jeveux de l'objet
     */
    BaseMesh( const std::string &name, const std::string &type );

    /**
     * @brief Read a Mesh file
     * @return return true if it succeeds, false otherwise
     */
    bool readMeshFile( const std::string &fileName, const std::string &format );

  public:
    /**
     * @typedef BaseMeshPtr
     * @brief Pointeur intelligent vers un BaseMesh
     */
    typedef std::shared_ptr< BaseMesh > BaseMeshPtr;

    /**
     * @brief Get the connectivity
     */
    const ConnectivityMeshExplorer &getConnectivityExplorer() const {
        _cellsType->updateValuePointer();
        _connectivity->build();
        return _explorer;
    };

    /**
     * @brief Return the connectivity
     */
    const JeveuxCollectionLong getConnectivity() const { return _connectivity; }

    const JeveuxCollectionLong getInverseConnectivity() const;

    /**
     * @brief Return the connectivity with MED numberings
     */
    const JeveuxCollectionLong getMedConnectivity() const;

    /**
     * @brief Return the MED type for each cell
     */
    const JeveuxVectorLong getMedCellsTypes() const;

    /**
     * @brief Recuperation des coordonnees du maillage
     * @return champ aux noeuds contenant les coordonnees des noeuds du maillage
     */
    MeshCoordinatesFieldPtr getCoordinates() const {
        _coordinates->updateValuePointers();
        return _coordinates;
    };

    ConstantFieldOnCellsRealPtr getCurvilinearAbscissa() const { return _curvAbsc; }

    /**
     * @brief Get all the names of group of cells
     * @return NamesMapChar24 _nameOfGrpCells
     */
    const NamesMapChar24 &getGroupsOfNodesMap() const { return _nameOfGrpCells; };

    /**
     * @brief Returns the number of nodes
     */
    ASTERINTEGER getNumberOfNodes() const;

    /**
     * @brief Returns the number of cells
     */
    ASTERINTEGER getNumberOfCells() const;

    /**
     * @brief Get all the names of nodes
     * @return NamesMapChar8 _nameOfNodes
     */
    const NamesMapChar8 &getNameOfNodesMap() const { return _nameOfNodes; };

    const NamesMapChar8 &getCellNameMap() const { return _nameOfCells; };

    std::string getNodeName( const ASTERINTEGER &index ) const;

    std::string getCellName( const ASTERINTEGER &index ) const;

    ASTERINTEGER getCellType( const ASTERINTEGER &index ) const;

    std::string getCellTypeName( const ASTERINTEGER &index ) const;

    /**
     * @brief Recuperation de la dimension du maillage
     */
    ASTERINTEGER getDimension() const;

    /**
     * @brief Teste l'existence d'un groupe de mailles dans le maillage
     * @return true si le groupe existe
     */
    virtual bool hasGroupOfCells( const std::string &name, const bool local = false ) const {
        AS_ASSERT( false );
        return false;
    };

    /**
     * @brief Teste l'existence d'un groupe de noeuds dans le maillage
     * @return true si le groupe existe
     */
    virtual bool hasGroupOfNodes( const std::string &name, const bool local = false ) const {
        AS_ASSERT( false );
        return false;
    };

    /**
     * @brief Returns the names of the groups of cells
     * @return VectorString
     */
    virtual VectorString getGroupsOfCells( const bool local = false ) const {
        AS_ASSERT( false );
        return {};
    };

    /**
     * @brief Returns the names of the groups of nodes
     * @return VectorString
     */
    virtual VectorString getGroupsOfNodes( const bool local = false ) const {
        AS_ASSERT( false );
        return {};
    };

    /**
     * @brief Returns the cells indexes of a group of cells
     * @return VectorLong
     */
    virtual VectorLong getCells( const std::string name ) const {
        AS_ASSERT( false );
        return {};
    }

    /**
     * @brief Returns the nodes indexes of a group of nodes
     * @return VectorLong
     */
    virtual VectorLong getNodes( const std::string name = std::string(),
                                 const bool localNumbering = true,
                                 const ASTERINTEGER same_rank = PythonBool::None ) const {
        AS_ASSERT( false );
        return {};
    }

    /**
     * @brief Returns the nodes indexes of a group of cells
     * @return VectorLong
     */
    virtual VectorLong getNodesFromCells( const std::string name, const bool localNumbering = true,
                                          const ASTERINTEGER same_rank = PythonBool::None ) const {
        AS_ASSERT( false );
        return {};
    }

    virtual VectorLong getNodesFromCells( const VectorLong &cells, const bool localNumbering = true,
                                          const ASTERINTEGER same_rank = PythonBool::None ) const {
        AS_ASSERT( false );
        return {};
    }

    /**
     * @brief Returns the nodes indexes of inner nodes
     * @return VectorLong
     */
    virtual VectorLong getInnerNodes() const {
        AS_ASSERT( false );
        return {};
    }

    /**
     * @brief Returns the nodes indexes of outer nodes
     * @return VectorLong
     */
    virtual VectorLong getOuterNodes() const {
        AS_ASSERT( false );
        return {};
    }

    /**
     * @brief Get the JeveuxVector for outer subdomain nodes
     * @return VectorLong
     */
    virtual const JeveuxVectorLong getNodesRank() const {
        AS_ASSERT( false );
        return {};
    }

    /**
     * @brief Get the JeveuxVector for outer subdomain cells
     * @return VectorLong
     */
    virtual const JeveuxVectorLong getCellsRank() const {
        AS_ASSERT( false );
        return {};
    }

    /**
     * @brief Fonction permettant de savoir si un maillage est vide (non relu par exemple)
     * @return retourne true si le maillage est vide
     */
    bool isEmpty() const { return !_dimensionInformations->exists(); };

    /**
     * @brief Fonction permettant de savoir si un maillage est parallel
     * @return retourne true si le maillage est parallel
     */
    virtual bool isParallel() const { return false; };

    /**
     * @brief Fonction permettant de savoir si un maillage est partiel
     * @return retourne true si le maillage est partiel
     */
    virtual bool isConnection() const { return false; };

    /**
     * @brief Tester le maillage a des cells quadratiques
     * @return true si quadratique
     */
    virtual bool isQuadratic() const {
        AS_ASSERT( false );
        return false;
    };

    /**
     * @brief Read a MED Mesh file
     * @return retourne true si tout est ok
     */
    virtual bool readMedFile( const std::string &fileName );

    /**
     * @brief Impression du maillage au format MED
     * @param fileName Nom du fichier MED ?? imprimer
     * @return true
     */
    bool printMedFile( const std::string fileName, bool local = true ) const;

    /**
     * @brief Get the mapping between local and global numbering of nodes
     * @return JeveuxVector of the indirection
     */
    virtual const JeveuxVectorLong getLocalToGlobalMapping() const {
        AS_ASSERT( false );
        return {};
    }

    /**
     * @brief Build the mesh
     * @return true if success
     */
    bool build();
};

/**
 * @typedef BaseMeshPtr
 * @brief Pointeur intelligent vers un BaseMesh
 */
typedef std::shared_ptr< BaseMesh > BaseMeshPtr;

#endif /* BASEMESH_H_ */
