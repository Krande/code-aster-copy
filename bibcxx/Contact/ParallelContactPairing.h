/**
 * @file ParallelContactPairing.h
 * @brief Fichier entete de la class ContactPairing
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

#ifdef ASTER_HAVE_MPI

#include "Contact/ContactPairing.h"
#include "Contact/ParallelContactNew.h"

class ParallelContactPairing : public ContactPairing {
  protected:
    /** @brief FiniteElementDescriptor on initial ParallelMeshPtr */
    ParallelContactFEDescriptorPtr _PFEDesc;
    FiniteElementDescriptorPtr _SFEDesc;
    ParallelContactNewPtr _pContNew;

  public:
    /**
     * @typedef ContactPairingPtr
     * @brief Pointeur intelligent vers un ContactPairing
     */
    typedef std::shared_ptr< ParallelContactPairing > ParallelContactPairingPtr;

    /**
     * @brief Constructeur
     */
    ParallelContactPairing() = delete;

    /**
     * @brief Constructeur
     */
    ParallelContactPairing( const std::string name, const ParallelContactNewPtr cont );

    /**
     * @brief Constructeur
     */
    ParallelContactPairing( const ParallelContactNewPtr cont )
        : ParallelContactPairing( ResultNaming::getNewResultName(), cont ) {};

    void buildFiniteElementDescriptor();

    ParallelContactFEDescriptorPtr getParallelFiniteElementDescriptor() const { return _PFEDesc; };
};

/**
 * @typedef ParallelContactPairingPtr
 * @brief Pointeur intelligent vers un ParallelContactPairing
 */
using ParallelContactPairingPtr = std::shared_ptr< ParallelContactPairing >;

#endif /* ASTER_HAVE_MPI */
