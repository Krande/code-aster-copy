
#include "astercxx.h"

#ifdef ASTER_HAVE_MPI

#ifndef PARALLELDOFNUMBERING_H_
#define PARALLELDOFNUMBERING_H_

/**
 * @file ParallelDOFNumbering.h
 * @brief Fichier entete de la classe ParallelDOFNumbering
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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

#include "Numbering/BaseDOFNumbering.h"
#include "Numbering/ParallelGlobalEquationNumbering.h"

#include <unordered_map>

/**
 * @class ParallelDOFNumbering
 * @brief Class definissant un nume_ddl_p
 * @author Nicolas Sellenet
 */
class ParallelDOFNumbering : public BaseDOFNumbering {
  private:
    /** @brief Objet '.NUME' */
    ParallelGlobalEquationNumberingPtr _globalNumbering;

    std::unordered_map< ASTERINTEGER, ASTERINTEGER > _global2localMap;

    /**
     * @brief Build the mapping from global to local numbering of the dof
     */
    void _buildGlobal2LocalMap() {
        getLocalToGlobalMapping()->updateValuePointer();
        ASTERINTEGER nloc = getLocalToGlobalMapping()->size();

        _global2localMap.reserve( nloc );
        for ( auto j = 0; j < nloc; j++ )
            _global2localMap[( *getLocalToGlobalMapping() )[j]] = j;
    };

  public:
    /**
     * @typedef ParallelDOFNumberingPtr
     * @brief Pointeur intelligent vers un ParallelDOFNumbering
     */
    typedef std::shared_ptr< ParallelDOFNumbering > ParallelDOFNumberingPtr;

    /**
     * @brief Constructeur
     */
    ParallelDOFNumbering();

    ParallelDOFNumbering( const std::string name, const FieldOnNodesDescriptionPtr fdof,
                          const ModelPtr model );

    /**
     * @brief Constructeur
     * @param name nom souhaité de la sd (utile pour le BaseDOFNumbering d'une sd_resu)
     */
    ParallelDOFNumbering( const std::string &name );

    /**
     * @brief Returns the GlobalEquationNumberingPtr
     */
    GlobalEquationNumberingPtr getGlobalNumbering() const { return _globalNumbering; };

    std::string getPhysicalQuantity() const { return _globalNumbering->getPhysicalQuantity(); };

    /**
     * @brief Methode permettant de savoir si l'objet est parallel
     * @return true
     */
    bool isParallel() { return true; };

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
     * @brief Get Rows Associated to all Ghost Dof
     */
    VectorLong getGhostRows( const bool local = false ) const;

    /**
     * @brief Get Rows owned locally (aka not Ghost)
     */
    VectorLong getNoGhostRows() const;

    /**
     * @brief Get Rows Associated to Lagrange Multipliers Dof
     */
    VectorLong getRowsAssociatedToLagrangeMultipliers( const bool local = false ) const;

    /**
     * @brief Get Assigned Components
     */
    VectorString getComponents() const;

    /**
     * @brief Get the mapping between local ang global numbering of the Dof
     */
    const JeveuxVectorLong getLocalToGlobalMapping() const;

    /**
     * @brief Return the local number of a global Dof
     * @return Return the local number if the row if present on the subdomain ; otherwise
     * raise an exception
     */
    const ASTERINTEGER globalToLocalRow( const ASTERINTEGER ) const;

    /**
     * @brief Return the global number of a local Dof
     * @return Return the global number if the row if present on the subdomain ; otherwise
     * raise an exception
     */
    const ASTERINTEGER localToGlobalRow( const ASTERINTEGER );

    /**
     * @brief Get model
     */
    ModelPtr getModel() const { return _globalNumbering->getModel(); };

    /**
     * @brief Set model
     */
    void setModel( const ModelPtr &model ) { _globalNumbering->setModel( model ); };

    /**
     * @brief Get mesh
     */
    BaseMeshPtr getMesh() const { return _globalNumbering->getMesh(); };
};

/**
 * @typedef ParallelDOFNumberingPtr
 * @brief Pointeur intelligent vers un ParallelDOFNumbering
 */
typedef std::shared_ptr< ParallelDOFNumbering > ParallelDOFNumberingPtr;

#endif /* PARALLELDOFNUMBERING_H_ */

#endif /* ASTER_HAVE_MPI */
