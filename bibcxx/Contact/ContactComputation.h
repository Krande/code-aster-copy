
/**
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

#pragma once

#include "astercxx.h"

#include "Contact/ContactNew.h"
#include "Contact/ContactPairing.h"
#include "DataFields/FieldOnCells.h"
#include "LinearAlgebra/ElementaryMatrix.h"

class ContactComputation {
  private:
    // Contact Definition
    ContactNewPtr _contact;

    /**
     * @brief Convert ELNO -> NOEU for virtual nodes
     */
    FieldOnNodesRealPtr convertVirtualField( const FieldOnCellsRealPtr field ) const;

  public:
    ContactComputation( const ContactNewPtr contact ) : _contact( contact ) {};

    /**
     * @brief Compute geometric gap
     */
    std::pair< FieldOnNodesRealPtr, FieldOnNodesRealPtr >
    geometricGap( const ContactPairingPtr pairing ) const;

    /**
     * @brief Compute contact mortar matrix
     */
    ElementaryMatrixDisplacementRealPtr contactMortarMatrix() const;

    FieldOnCellsRealPtr contactData( const ContactPairingPtr pairing, const MaterialFieldPtr mater,
                                     const bool &initial_contact ) const;

    /**
     * @brief Compute contact coefficient field (COEF_CONT)
     */
    std::pair< FieldOnNodesRealPtr, FieldOnNodesRealPtr > contactCoefficient() const;
};

/**
 * @typedef ContactComputationPtr
 * @brief Pointeur intelligent vers un ContactComputation
 */
typedef std::shared_ptr< ContactComputation > ContactComputationPtr;
