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

#ifndef DOFNUMBERING_H_
#define DOFNUMBERING_H_

/**
 * @file DOFNumbering.h
 * @brief Fichier entete de la classe DOFNumbering
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "astercxx.h"

#include "LinearAlgebra/ElementaryMatrix.h"
#include "LinearAlgebra/MatrixStorage.h"
#include "Loads/DirichletBC.h"
#include "Loads/ListOfLoads.h"
#include "Loads/MechanicalLoad.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/BaseMesh.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Numbering/BaseDOFNumbering.h"
#include "Numbering/FieldOnNodesDescription.h"

/**
 * @class DOFNumbering
 * @brief Class definissant un nume_ddl
 * @author Nicolas Sellenet
 */
class DOFNumbering : public BaseDOFNumbering {
  private:
    /** @brief Objet '.NUME' */
    GlobalEquationNumberingPtr _globalNumbering;
    /** @brief Objet '.NUML' */
    LocalEquationNumberingPtr _localNumbering;

  public:
    /**
     * @typedef DOFNumberingPtr
     * @brief Pointeur intelligent vers un DOFNumbering
     */
    typedef std::shared_ptr< DOFNumbering > DOFNumberingPtr;

    /**
     * @brief Constructeur
     */
    DOFNumbering();

    DOFNumbering( const std::string name, const FieldOnNodesDescriptionPtr fdof,
                  const ModelPtr model );

    /**
     * @brief Constructeur
     * @param name nom souhaité de la sd (utile pour le BaseDOFNumbering d'une sd_resu)
     */
    DOFNumbering( const std::string name );

    /**
     * @brief Returns the GlobalEquationNumberingPtr
     */
    virtual GlobalEquationNumberingPtr getGlobalNumbering() const { return _globalNumbering; };

    /**
     * @brief Get Physical Quantity
     */
    std::string getPhysicalQuantity() const;

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
     * @brief Get The Components Associated To A Given Node
     */
    VectorString getComponentsAssociatedToNode( const ASTERINTEGER node,
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
     * @brief get the Row index Associated To the Component of a Node
     */
    ASTERINTEGER getRowAssociatedToNodeComponent( const ASTERINTEGER node, const std::string comp,
                                                  const bool local = false ) const;

    /**
     * @brief Get Rows Associated to all Physical Dof
     */
    VectorLong getRowsAssociatedToPhysicalDofs( const bool local = false ) const;

    /**
     * @brief Get Rows Associated to Lagrange Multipliers Dof
     */
    VectorLong getRowsAssociatedToLagrangeMultipliers( const bool local = false ) const;

    /**
     * @brief Get Assigned Components
     */
    VectorString getComponents() const;

    /**
     * @brief Get model
     */
    ModelPtr getModel() const { return _globalNumbering->getModel(); };

    /**
     * @brief Set model
     */
    bool setModel( const ModelPtr &model ) { return _globalNumbering->setModel( model ); };
};

/**
 * @typedef DOFNumberingPtr
 * @brief Enveloppe d'un pointeur intelligent vers un DOFNumbering
 * @author Nicolas Sellenet
 */
typedef std::shared_ptr< DOFNumbering > DOFNumberingPtr;

#endif /* DOFNUMBERING_H_ */
