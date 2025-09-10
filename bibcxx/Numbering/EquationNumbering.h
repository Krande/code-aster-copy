/**
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
#pragma once

#include "astercxx.h"

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/BaseMesh.h"
#include "Modeling/Model.h"

class BaseEquationNumbering : public DataStructure {

  public:
    BaseEquationNumbering( const std::string baseName = DataStructureNaming::getNewName() )
        : DataStructure( baseName, 19, "NUME_EQUA" ) {};

    /**
     * @brief Returns the vector of local to global numbering
     */
    virtual const JeveuxVectorLong getLocalToGlobal() const {
        throw std::runtime_error( "Vector LocalToGlobal doesn't exist in sequential" );
        return JeveuxVectorLong( "RIEN" );
    };

    /**
     * @brief Returns the vector of the rank owning the local dof number
     */
    virtual const JeveuxVectorLong getLocalToRank() const {
        throw std::runtime_error( "Vector LocalToRank doesn't exist in sequential" );
        return JeveuxVectorLong( "RIEN" );
    };

    /**
     * @brief Are Lagrange Multipliers used for BC or MPC
     */
    virtual bool useLagrangeDOF() const = 0;

    /**
     * @brief Are Single Lagrange Multipliers used for BC or MPC
     */
    virtual bool useSingleLagrangeDOF() const = 0;

    /**
     * @brief Get The Component Associated To A Given Row
     */
    virtual std::string getComponentFromDOF( const ASTERINTEGER dof,
                                             const bool local = false ) const = 0;
    /**
     * @brief Get The Node Id Associated To A Given Row
     */
    virtual ASTERINTEGER getNodeFromDOF( const ASTERINTEGER dof,
                                         const bool local = false ) const = 0;

    /**
     * @brief Return true if a physical dof is Associated To A Given Row
     */
    virtual bool isPhysicalDOF( const ASTERINTEGER dof, const bool local = false ) const = 0;

    /**
     * @brief Get The total number of Dofs
     */
    virtual ASTERINTEGER getNumberOfDOFs( const bool local = false ) const = 0;

    /**
     * @brief Get Rows Associated to all Physical Dof
     */
    virtual VectorLong getPhysicalDOFs( const bool local = false ) const = 0;

    /**
     * @brief Get Rows Associated to Lagrange Multipliers Dof
     */
    virtual VectorLong getLagrangeDOFs( const bool local = false ) const = 0;

    virtual ASTERINTEGER getDOFFromNodeAndComponent( const ASTERINTEGER &node,
                                                     const std::string &comp,
                                                     const bool local = false ) const = 0;

    /**
     * @brief Get Rows Associated to all Ghost Dof
     */
    VectorLong getGhostDOFs( const bool local = false, const bool lastLayerOnly = false ) const {
        throw std::runtime_error( "No ghost DOF in sequential" );
        return VectorLong();
    };

    /**
     * @brief Get Rows owned locally (aka not Ghost)
     */
    VectorLong getNoGhostDOFs() const {
        throw std::runtime_error( "No ghost DOF in sequential" );
        return VectorLong();
    };

    /**
     * @brief Return the local number of a global Dof
     * @return Return the local number if the dof if present on the subdomain ; otherwise
     * raise an exception
     */
    ASTERINTEGER globalToLocalDOF( const ASTERINTEGER ) const {
        throw std::runtime_error( "Vector globalToLocalDOF doesn't exist in sequential" );
        return -1;
    };

    /**
     * @brief Return the global number of a local Dof
     * @return Return the global number if the dof if present on the subdomain ; otherwise
     * raise an exception
     */
    ASTERINTEGER localToGlobalDOF( const ASTERINTEGER ) const {
        throw std::runtime_error( "Vector globalToLocalDOF doesn't exist in sequential" );
        return -1;
    };
};

/**
 * @class EquationNumbering
 * @brief Class definissant un NUME_EQUA
 */
class EquationNumbering : public BaseEquationNumbering {
  protected:
    /** @brief Objet Jeveux '.NEQU' */
    JeveuxVectorLong _numberOfEquations;
    /** @brief Objet Jeveux '.REFN' */
    JeveuxVectorChar24 _informations;
    /** @brief Objet Jeveux '.DELG' */
    JeveuxVectorLong _lagrangianInformations;
    /** @brief Objet Jeveux '.PRNO' */
    JeveuxCollectionLong _componentsOnNodes;
    /** @brief Objet Jeveux '.LILI' */
    NamesMapChar24 _namesOfGroupOfCells;
    /** @brief Objet Jeveux '.NUEQ' */
    JeveuxVectorLong _indexationVector;
    /** @brief Objet Jeveux '.DEEQ' */
    JeveuxVectorLong _nodeAndComponentsIdFromDOF;
    /** @brief Mesh */
    BaseMeshPtr _mesh;
    /** @brief Model */
    ModelPtr _model;

    std::map< ASTERINTEGER, std::string > _componentsNumber2Name;

    /**
     * @brief Build the mapping between the component number to its name
     */
    void _buildAllComponentsId2Name();

  public:
    /**
     * @typedef EquationNumberingPtr
     * @brief Pointeur intelligent vers un EquationNumbering
     */
    typedef std::shared_ptr< EquationNumbering > EquationNumberingPtr;

    EquationNumbering( const std::string &baseName );

    EquationNumbering();

    /**
     * @brief Surcharge de l'operateur =
     */
    bool operator==( EquationNumbering &toCompare );

    bool operator!=( EquationNumbering &toCompare ) { return !( *this == toCompare ); }

    /**
     * @brief Returns a vector of information of the Lagrange multipliers
     */
    const JeveuxVectorLong getLagrangianInformations() const { return _lagrangianInformations; }

    /**
     * @brief Returns a vector of information on the numer of equations
     */
    const JeveuxVectorLong getNumberOfEquations() const { return _numberOfEquations; }

    /**
     * @brief Get model
     */
    ModelPtr getModel() const { return _model; };

    /**
     * @brief Set model
     */
    void setModel( const ModelPtr &model );

    /**
     * @brief Get Mesh
     */
    BaseMeshPtr getMesh() const { return _mesh; };

    /**
     * @brief Set Mesh
     */
    void setMesh( const BaseMeshPtr &mesh );

    /**
     * @brief Get Physical Quantity
     */
    std::string getPhysicalQuantity() const;

    bool isParallel() const { return false; };

    bool exists() const;

    /**
     * @brief Returns a vector with node index for each DOFs
     */
    VectorLong getNodesFromDOF() const;

    /**
     * @brief Return list of DOFs
     * @param sameRank True: Use only owned nodes / False: Use all nodes
     * @param list_cmp empty: Use all cmp / keep only cmp given
     * @param groupsOfCells empty: Use all nodes / keep only nodes given
     */
    VectorLong getDOFs( const bool sameRank = false, const VectorString &list_cmp = {},
                        const VectorString &list_grpno = {} ) const;

    /**
     * @brief Returns a vector with node index and component name for each DOFs
     */
    VectorPairLong getNodeAndComponentIdFromDOF( const bool local = true ) const;

    PairLong getNodeAndComponentIdFromDOF( const ASTERINTEGER dof, const bool local = true ) const;

    /**
     * @brief Returns a vector with node index and component name for each DOFs
     */
    std::vector< std::pair< ASTERINTEGER, std::string > >
    getNodeAndComponentFromDOF( const bool local = true ) const;

    /**
     * @brief Returns a pair with node index and component name for given DOF
     */
    std::pair< ASTERINTEGER, std::string >
    getNodeAndComponentFromDOF( const ASTERINTEGER dof, const bool local = true ) const;

    /**
     * @brief Maps between node id and name of components to DOF
     */
    std::map< PairLong, ASTERINTEGER >
    getDOFFromNodeAndComponentId( const bool local = true ) const;

    std::map< std::pair< ASTERINTEGER, std::string >, ASTERINTEGER >
    getDOFFromNodeAndComponent( const bool local = true ) const;

    ASTERINTEGER getDOFFromNodeAndComponent( const ASTERINTEGER &node, const std::string &comp,
                                             const bool local = false ) const;
    /**
     * @brief Get componants
     */
    VectorString getComponents() const;
    SetLong getComponentsId() const;

    /**
     * @brief Maps between name of components and the number
     */
    std::map< std::string, ASTERINTEGER > getComponentsNameToId() const;
    std::map< ASTERINTEGER, std::string > getComponentsIdToName() const;

    /**
     * @brief Are Lagrange Multipliers used for BC or MPC
     */
    bool useLagrangeDOF() const;

    /**
     * @brief Get Rows Associated to all Ghost Dof
     */
    VectorLong getGhostDOFs( const bool local = false, const bool lastLayerOnly = false ) const {
        // throw std::runtime_error( "No ghost DOF in sequential" );
        return VectorLong();
    };

    /**
     * @brief Get Rows owned locally (aka not Ghost)
     */
    VectorLong getNoGhostDOFs( const bool local = false ) const;

    /**
     * @brief Are Single Lagrange Multipliers used for BC or MPC
     */
    bool useSingleLagrangeDOF() const;

    /**
     * @brief Get The Component Associated To A Given Row
     */
    std::string getComponentFromDOF( const ASTERINTEGER dof, const bool local = false ) const;
    /**
     * @brief Get The Node Id Associated To A Given Row
     */
    ASTERINTEGER getNodeFromDOF( const ASTERINTEGER dof, const bool local = false ) const;

    /**
     * @brief Get The Rows Associated To A Given Node Id
     */
    VectorLong getDOFsFromNode( const ASTERINTEGER node, const bool local = false ) const;

    /**
     * @brief Return true if a physical dof is Associated To A Given Row
     */
    bool isPhysicalDOF( const ASTERINTEGER dof, const bool local = false ) const;

    /**
     * @brief Get The total number of Dofs
     */
    ASTERINTEGER getNumberOfDOFs( const bool local = false ) const;

    /**
     * @brief Get Rows Associated to all Physical Dof
     */
    VectorLong getPhysicalDOFs( const bool local = false ) const;

    /**
     * @brief Get Rows Associated to Lagrange Multipliers Dof
     */
    VectorLong getLagrangeDOFs( const bool local = false ) const;

    /**
     * @brief Get Rows Associated to first and second Lagrange Multipliers Dof
     */

    std::map< ASTERINTEGER, VectorLong > getDictOfLagrangeDOFs( const bool local = false ) const;

    /**
     * @brief Get (Nodes And Components) and dofs corresponding to component name and list of node
     * groups
     */
    std::pair< std::pair< VectorLong, VectorString >, VectorLong >
    getDOFsWithDescription( const VectorString &, const VectorString &, const bool local = false,
                            const ASTERINTEGER same_rank = PythonBool::None ) const;

    /**
     * @brief Get (Nodes And Components) and dofs corresponding to component name and list of node
     * groups
     */
    std::pair< std::pair< VectorLong, VectorString >, VectorLong >
    getDOFsWithDescription( const VectorString &, const VectorLong &, const bool local = false,
                            const ASTERINTEGER same_rank = PythonBool::None ) const;

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
     */
    void updateValuePointers();
};

/**
 * @typedef EquationNumberingPtr
 * @brief Enveloppe d'un pointeur intelligent vers un EquationNumbering
 */
typedef std::shared_ptr< EquationNumbering > EquationNumberingPtr;
