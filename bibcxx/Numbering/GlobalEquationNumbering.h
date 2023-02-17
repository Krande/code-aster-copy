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

/**
 * @class GlobalEquationNumbering
 * @brief Class definissant un NUME_EQUA
 */
class GlobalEquationNumbering : public DataStructure {
  private:
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

  public:
    /**
     * @typedef GlobalEquationNumberingPtr
     * @brief Pointeur intelligent vers un GlobalEquationNumbering
     */
    typedef std::shared_ptr< GlobalEquationNumbering > GlobalEquationNumberingPtr;

    GlobalEquationNumbering( const std::string &baseName );

    GlobalEquationNumbering();

    /**
     * @brief Surcharge de l'operateur =
     */
    bool operator==( GlobalEquationNumbering &toCompare );

    bool operator!=( GlobalEquationNumbering &toCompare ) { return !( *this == toCompare ); }

    /**
     * @brief Returns a vector of information of the Lagrange multipliers
     */
    const JeveuxVectorLong getLagrangianInformations() const { return _lagrangianInformations; }

    /**
     * @brief Returns a vector of information on the numer of equations
     */
    const JeveuxVectorLong getNumberOfEquations() const { return _numberOfEquations; }

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
     * @brief Returns number of DOFs
     */
    ASTERINTEGER getNumberOfDofs() const;

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
    SetString getComponents() const;
    SetLong getComponentsNumber() const;

    /**
     * @brief Maps between name of components and the number
     */
    std::map< std::string, ASTERINTEGER > getComponentsName2Number() const;
    std::map< ASTERINTEGER, std::string > getComponentsNumber2Name() const;

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
     */
    void updateValuePointers();
};

/**
 * @typedef GlobalEquationNumberingPtr
 * @brief Enveloppe d'un pointeur intelligent vers un GlobalEquationNumbering
 */
typedef std::shared_ptr< GlobalEquationNumbering > GlobalEquationNumberingPtr;
