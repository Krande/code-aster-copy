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
#include "Numbering/GlobalEquationNumbering.h"

#pragma once

/**
 * @class ParallelGlobalEquationNumbering
 * @brief Class definissant un NUME_EQUA
 */
class ParallelGlobalEquationNumbering : public GlobalEquationNumbering {
  protected:
    /** @brief Objet Jeveux '.NULG' */
    JeveuxVectorLong _localToGlobal;
    /** @brief Objet Jeveux '.PDDL' */
    JeveuxVectorLong _localToRank;

    std::unordered_map< ASTERINTEGER, ASTERINTEGER > _global2localMap;

    /**
     * @brief Build the mapping from global to local numbering of the dof
     */
    void _buildGlobal2LocalMap();

  public:
    /**
     * @typedef GlobalEquationNumberingPtr
     * @brief Pointeur intelligent vers un ParallelGlobalEquationNumbering
     */
    typedef std::shared_ptr< ParallelGlobalEquationNumbering > ParallelGlobalEquationNumberingPtr;

    ParallelGlobalEquationNumbering( const std::string &baseName );

    /**
     * @brief Returns the vector of local to global numbering
     */
    const JeveuxVectorLong getLocalToGlobal() const { return _localToGlobal; };

    /**
     * @brief Returns the vector of the rank owning the local dof number
     */
    const JeveuxVectorLong getLocalToRank() const { return _localToRank; };

    /**
     * @brief Returns a vector with node index and component name for each DOFs
     */
    VectorPairLong getNodesAndComponentsNumberFromDOF( const bool local = true ) const;

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
    const JeveuxVectorLong getLocalToGlobalMapping() const { return getLocalToGlobal(); };

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
     * @brief Are Lagrange Multipliers used for BC or MPC
     */
    bool useLagrangeMultipliers() const;

    /**
     * @brief Are Single Lagrange Multipliers used for BC or MPC
     */
    bool useSingleLagrangeMultipliers() const;

    bool isParallel() const { return true; };
};

/**
 * @typedef ParallelGlobalEquationNumberingPtr
 * @brief Enveloppe d'un pointeur intelligent vers un ParallelGlobalEquationNumbering
 */
typedef std::shared_ptr< ParallelGlobalEquationNumbering > ParallelGlobalEquationNumberingPtr;
