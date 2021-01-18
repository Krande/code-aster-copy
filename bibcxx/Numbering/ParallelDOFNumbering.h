
#include "astercxx.h"

#ifdef ASTER_HAVE_MPI

#ifndef PARALLELDOFNUMBERING_H_
#define PARALLELDOFNUMBERING_H_

/**
 * @file ParallelDOFNumbering.h
 * @brief Fichier entete de la classe ParallelDOFNumbering
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "Numbering/DOFNumbering.h"

/**
 * @class ParallelDOFNumberingClass
 * @brief Cette classe decrit un maillage Aster parallèle
 * @author Nicolas Sellenet
 */
class ParallelDOFNumberingClass : public BaseDOFNumberingClass {
  private:
  public:
    /**
     * @typedef ParallelDOFNumberingPtr
     * @brief Pointeur intelligent vers un ParallelDOFNumberingClass
     */
    typedef boost::shared_ptr< ParallelDOFNumberingClass > ParallelDOFNumberingPtr;

    /**
     * @brief Constructeur
     */
    ParallelDOFNumberingClass( const JeveuxMemory memType = Permanent )
        : BaseDOFNumberingClass( "NUME_DDL_P", memType ){};

    /**
     * @brief Constructeur
     * @param name nom souhaité de la sd (utile pour le BaseDOFNumberingClass d'une sd_resu)
     */
    ParallelDOFNumberingClass( const std::string name, const JeveuxMemory memType = Permanent )
        : BaseDOFNumberingClass( name, "NUME_DDL_P", memType ){};

    /**
     * @brief Methode permettant de savoir si l'objet est parallel
     * @return true
     */
    bool isParallel() { return true; };

    /**
     * @brief Methode permettant de definir les matrices elementaires
     * @param currentMatrix objet ElementaryMatrix
     */
    void setElementaryMatrix( const ElementaryMatrixDisplacementRealPtr &currentMatrix )

    {
        if ( !currentMatrix->getModel()->getMesh()->isParallel() )
            throw std::runtime_error( "Mesh must be parallel" );
        BaseDOFNumberingClass::setElementaryMatrix( currentMatrix );
    };

    /**
     * @brief Methode permettant de definir les matrices elementaires
     * @param currentMatrix objet ElementaryMatrix
     */
    void setElementaryMatrix( const ElementaryMatrixDisplacementComplexPtr &currentMatrix )

    {
        if ( !currentMatrix->getModel()->getMesh()->isParallel() )
            throw std::runtime_error( "Mesh must be parallel" );
        BaseDOFNumberingClass::setElementaryMatrix( currentMatrix );
    };

    /**
     * @brief Methode permettant de definir les matrices elementaires
     * @param currentMatrix objet ElementaryMatrix
     */
    void setElementaryMatrix( const ElementaryMatrixTemperatureRealPtr &currentMatrix )

    {
        if ( !currentMatrix->getModel()->getMesh()->isParallel() )
            throw std::runtime_error( "Mesh must be parallel" );
        BaseDOFNumberingClass::setElementaryMatrix( currentMatrix );
    };

    /**
     * @brief Methode permettant de definir les matrices elementaires
     * @param currentMatrix objet ElementaryMatrix
     */
    void setElementaryMatrix( const ElementaryMatrixPressureComplexPtr &currentMatrix )

    {
        if ( !currentMatrix->getModel()->getMesh()->isParallel() )
            throw std::runtime_error( "Mesh must be parallel" );
        BaseDOFNumberingClass::setElementaryMatrix( currentMatrix );
    };

    /**
     * @brief Methode permettant de definir le modele
     * @param currentModel Modele de la numerotation
     */
    void setModel( const ModelPtr &currentModel ) {
        if ( !currentModel->getMesh()->isParallel() )
            throw std::runtime_error( "Mesh must be parallel" );
        BaseDOFNumberingClass::setModel( currentModel );
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
    ASTERINTEGER getRowAssociatedToNodeComponent(const ASTERINTEGER node, const std::string comp)
                                                                                            const {
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
 * @typedef ParallelDOFNumberingPtr
 * @brief Pointeur intelligent vers un ParallelDOFNumberingClass
 */
typedef boost::shared_ptr< ParallelDOFNumberingClass > ParallelDOFNumberingPtr;

#endif /* PARALLELDOFNUMBERING_H_ */

#endif /* ASTER_HAVE_MPI */
