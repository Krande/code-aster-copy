/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org             */
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
 *   Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
#include "boost/variant.hpp"

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


/**
 * @class DOFNumbering
 * @brief Class definissant un nume_ddl
 * @author Nicolas Sellenet
 */
class DOFNumbering : public BaseDOFNumbering {
  public:
    /**
     * @typedef DOFNumberingPtr
     * @brief Pointeur intelligent vers un DOFNumbering
     */
    typedef boost::shared_ptr< DOFNumbering > DOFNumberingPtr;

    /**
     * @brief Constructeur
     */
    DOFNumbering( const JeveuxMemory memType = Permanent )
        : BaseDOFNumbering( "NUME_DDL", memType ){};

    /**
     * @brief Constructeur
     * @param name nom souhaité de la sd (utile pour le BaseDOFNumbering d'une sd_resu)
     */
    DOFNumbering( const std::string name, const JeveuxMemory memType = Permanent )
        : BaseDOFNumbering( name, "NUME_DDL", memType ){};

    /**
     * @brief Methode permettant de definir les matrices elementaires
     * @param currentMatrix objet ElementaryMatrix
     */
    void setElementaryMatrix( const ElementaryMatrixDisplacementRealPtr &currentMatrix )

    {
        if ( currentMatrix->getModel()->getMesh()->isParallel() )
            throw std::runtime_error( "Mesh must not be parallel" );
        BaseDOFNumbering::setElementaryMatrix( currentMatrix );
    };

    /**
     * @brief Methode permettant de definir les matrices elementaires
     * @param currentMatrix objet ElementaryMatrix
     */
    void setElementaryMatrix( const ElementaryMatrixDisplacementComplexPtr &currentMatrix )

    {
        if ( currentMatrix->getModel()->getMesh()->isParallel() )
            throw std::runtime_error( "Mesh must not be parallel" );
        BaseDOFNumbering::setElementaryMatrix( currentMatrix );
    };

    /**
     * @brief Methode permettant de definir les matrices elementaires
     * @param currentMatrix objet ElementaryMatrix
     */
    void setElementaryMatrix( const ElementaryMatrixTemperatureRealPtr &currentMatrix )

    {
        if ( currentMatrix->getModel()->getMesh()->isParallel() )
            throw std::runtime_error( "Mesh must not be parallel" );
        BaseDOFNumbering::setElementaryMatrix( currentMatrix );
    };

    /**
     * @brief Methode permettant de definir les matrices elementaires
     * @param currentMatrix objet ElementaryMatrix
     */
    void setElementaryMatrix( const ElementaryMatrixPressureComplexPtr &currentMatrix )

    {
        if ( currentMatrix->getModel()->getMesh()->isParallel() )
            throw std::runtime_error( "Mesh must not be parallel" );
        BaseDOFNumbering::setElementaryMatrix( currentMatrix );
    };

    /**
     * @brief Methode permettant de definir le modele
     * @param currentModel Modele de la numerotation
     */
    void setModel( const ModelPtr &currentModel ) {
        if ( currentModel->getMesh()->isParallel() )
            throw std::runtime_error( "Mesh must not be parallel" );
        BaseDOFNumbering::setModel( currentModel );
    };

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
    std::string getComponentAssociatedToRow(const ASTERINTEGER row, const bool local) const;
    std::string getComponentAssociatedToRow(const ASTERINTEGER row) const {
        return getComponentAssociatedToRow(row, false);
    };

    /**
     * @brief Get The Components Associated To A Given Node
     */
    VectorString getComponentsAssociatedToNode(const ASTERINTEGER node, const bool local) const;
    VectorString getComponentsAssociatedToNode(const ASTERINTEGER node) const {
        return getComponentsAssociatedToNode(node, false);
    };

    /**
     * @brief Get The Node Id Associated To A Given Row
     */
    ASTERINTEGER getNodeAssociatedToRow(const ASTERINTEGER row, const bool local) const;
    ASTERINTEGER getNodeAssociatedToRow(const ASTERINTEGER row) const {
        return getNodeAssociatedToRow(row, false);
    };

    /**
     * @brief Get The total number of Dofs
     */
    ASTERINTEGER getNumberOfDofs(const bool local) const;
    ASTERINTEGER getNumberOfDofs() const {
        return getNumberOfDofs(false);
    };

    /**
     * @brief get the Row index Associated To the Component of a Node
     */
    ASTERINTEGER getRowAssociatedToNodeComponent(const ASTERINTEGER node, const std::string comp,
                                                                          const bool local) const;
    ASTERINTEGER getRowAssociatedToNodeComponent(const ASTERINTEGER node,
                                                 const std::string comp) const {
        return getRowAssociatedToNodeComponent(node, comp, false);
    };

    /**
     * @brief Get Rows Associated to all Physical Dof
     */
    VectorLong getRowsAssociatedToPhysicalDofs(const bool local) const;
    VectorLong getRowsAssociatedToPhysicalDofs() const {
        return getRowsAssociatedToPhysicalDofs(false);
    };

    /**
     * @brief Get Rows Associated to Lagrange Multipliers Dof
     */
    VectorLong getRowsAssociatedToLagrangeMultipliers(const bool local) const;
    VectorLong getRowsAssociatedToLagrangeMultipliers() const {
        return getRowsAssociatedToLagrangeMultipliers(false);
    };

    /**
     * @brief Get Assigned Components
     */
    VectorString getComponents() const;


};

/**
 * @typedef BaseDOFNumberingPtr
 * @brief Enveloppe d'un pointeur intelligent vers un BaseDOFNumbering
 * @author Nicolas Sellenet
 */
typedef boost::shared_ptr< BaseDOFNumbering > BaseDOFNumberingPtr;

/**
 * @typedef DOFNumberingPtr
 * @brief Enveloppe d'un pointeur intelligent vers un DOFNumbering
 * @author Nicolas Sellenet
 */
typedef boost::shared_ptr< DOFNumbering > DOFNumberingPtr;

#endif /* DOFNUMBERING_H_ */
