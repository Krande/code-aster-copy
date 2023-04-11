/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

/**
 * @section LICENCE
 *   Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
 *   This file is part of code_aster.
 *
 *   code_aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   code_aster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 */

#include "astercxx.h"

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/BaseMesh.h"
#include "Modeling/Model.h"

#pragma once

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
    virtual bool useLagrangeMultipliers() const = 0;

    /**
     * @brief Are Single Lagrange Multipliers used for BC or MPC
     */
    virtual bool useSingleLagrangeMultipliers() const = 0;

    /**
     * @brief Get The Component Associated To A Given Row
     */
    virtual std::string getComponentAssociatedToRow( const ASTERINTEGER row,
                                                     const bool local = false ) const = 0;
    /**
     * @brief Get The Node Id Associated To A Given Row
     */
    virtual ASTERINTEGER getNodeAssociatedToRow( const ASTERINTEGER row,
                                                 const bool local = false ) const = 0;

    /**
     * @brief Return true if a physical dof is Associated To A Given Row
     */
    virtual bool isRowAssociatedToPhysical( const ASTERINTEGER row,
                                            const bool local = false ) const = 0;

    /**
     * @brief Get The total number of Dofs
     */
    virtual ASTERINTEGER getNumberOfDofs( const bool local = false ) const = 0;

    /**
     * @brief Get Rows Associated to all Physical Dof
     */
    virtual VectorLong getRowsAssociatedToPhysicalDofs( const bool local = false ) const = 0;

    /**
     * @brief Get Rows Associated to Lagrange Multipliers Dof
     */
    virtual VectorLong getRowsAssociatedToLagrangeMultipliers( const bool local = false ) const = 0;

    /**
     * @brief Get Rows Associated to all Ghost Dof
     */
    VectorLong getGhostRows( const bool local = false ) const {
        throw std::runtime_error( "Vector LocagetGhostRowslToRank doesn't exist in sequential" );
        return VectorLong();
    };

    /**
     * @brief Get Rows owned locally (aka not Ghost)
     */
    VectorLong getNoGhostRows() const {
        throw std::runtime_error( "Vector LocagetGhostRowslToRank doesn't exist in sequential" );
        return VectorLong();
    };

    /**
     * @brief Return the local number of a global Dof
     * @return Return the local number if the row if present on the subdomain ; otherwise
     * raise an exception
     */
    const ASTERINTEGER globalToLocalRow( const ASTERINTEGER ) const {
        throw std::runtime_error( "Vector globalToLocalRow doesn't exist in sequential" );
        return -1;
    };

    /**
     * @brief Return the global number of a local Dof
     * @return Return the global number if the row if present on the subdomain ; otherwise
     * raise an exception
     */
    const ASTERINTEGER localToGlobalRow( const ASTERINTEGER ) {
        throw std::runtime_error( "Vector globalToLocalRow doesn't exist in sequential" );
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
    JeveuxVectorLong _nodeAndComponentsNumberFromDOF;
    /** @brief Mesh */
    BaseMeshPtr _mesh;
    /** @brief Model */
    ModelPtr _model;

    std::map< ASTERINTEGER, std::string > _componentsNumber2Name;

    /**
     * @brief Build the mapping between the component number to its name
     */
    void _buildAllComponentsNumber2Name();

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
                        const VectorLong &list_nodes = {} ) const;

    /**
     * @brief Returns a vector with node index and component name for each DOFs
     */
    VectorPairLong getNodesAndComponentsNumberFromDOF( const bool local = true ) const;

    PairLong getNodeAndComponentNumberFromDOF( const ASTERINTEGER dof,
                                               const bool local = true ) const;

    /**
     * @brief Returns a vector with node index and component name for each DOFs
     */
    std::vector< std::pair< ASTERINTEGER, std::string > >
    getNodesAndComponentsFromDOF( const bool local = true ) const;
    std::pair< ASTERINTEGER, std::string >
    getNodeAndComponentFromDOF( const ASTERINTEGER dof, const bool local = true ) const;

    /**
     * @brief Maps between node id and name of components to DOF
     */
    std::map< PairLong, ASTERINTEGER >
    getDOFsFromNodesAndComponentsNumber( const bool local = true ) const;

    std::map< std::pair< ASTERINTEGER, std::string >, ASTERINTEGER >
    getDOFsFromNodesAndComponents( const bool local = true ) const;

    /**
     * @brief Get componants
     */
    VectorString getComponents() const;
    SetLong getComponentsNumber() const;

    /**
     * @brief Maps between name of components and the number
     */
    std::map< std::string, ASTERINTEGER > getComponentsName2Number() const;
    std::map< ASTERINTEGER, std::string > getComponentsNumber2Name() const;

    /**
     * @brief Are Lagrange Multipliers used for BC or MPC
     */
    bool useLagrangeMultipliers() const;

    /**
     * @brief Are Single Lagrange Multipliers used for BC or MPC
     */
    bool useSingleLagrangeMultipliers() const;

    /**
     * @brief Get The Component Associated To A Given Row
     */
    std::string getComponentAssociatedToRow( const ASTERINTEGER row,
                                             const bool local = false ) const;
    /**
     * @brief Get The Node Id Associated To A Given Row
     */
    ASTERINTEGER getNodeAssociatedToRow( const ASTERINTEGER row, const bool local = false ) const;

    /**
     * @brief Return true if a physical dof is Associated To A Given Row
     */
    bool isRowAssociatedToPhysical( const ASTERINTEGER row, const bool local = false ) const;

    /**
     * @brief Get The total number of Dofs
     */
    ASTERINTEGER getNumberOfDofs( const bool local = false ) const;

    /**
     * @brief Get Rows Associated to all Physical Dof
     */
    VectorLong getRowsAssociatedToPhysicalDofs( const bool local = false ) const;

    /**
     * @brief Get Rows Associated to Lagrange Multipliers Dof
     */
    VectorLong getRowsAssociatedToLagrangeMultipliers( const bool local = false ) const;

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
