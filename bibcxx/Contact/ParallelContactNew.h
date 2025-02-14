/**
 * @file ParallelContactNew.h
 * @brief Fichier entete de la class ContactNew
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

#include "Contact/ContactNew.h"
#include "Modeling/ParallelContactFEDescriptor.h"

#ifdef ASTER_HAVE_MPI

class ParallelContactNew : public ContactNew {
  protected:
    /** @brief Mesh */
    ParallelMeshPtr _pMesh;
    /** @brief Parallel model */
    ModelPtr _pModel;
    /** @brief FiniteElementDescriptor on initial ParallelMeshPtr */
    ParallelContactFEDescriptorPtr _PFEDesc;
    FiniteElementDescriptorPtr _SFEDesc;

    /**
     * @brief Constructeur
     */
    ParallelContactNew( const std::string name, const ModelPtr model, ParallelMeshPtr mesh,
                        const std::string type );

  public:
    /**
     * @typedef ContactNewPtr
     * @brief Pointeur intelligent vers un ContactNew
     */
    typedef std::shared_ptr< ParallelContactNew > ParallelContactNewPtr;

    /**
     * @brief Constructeur
     */
    ParallelContactNew() = delete;

    /**
     * @brief Constructeur
     */
    ParallelContactNew( const std::string name, const ModelPtr model, ParallelMeshPtr mesh )
        : ParallelContactNew( name, model, mesh, "CHAR_CONT" ) {};

    /**
     * @brief Constructeur
     */
    ParallelContactNew( const ModelPtr model, ParallelMeshPtr mesh )
        : ParallelContactNew( ResultNaming::getNewResultName(), model, mesh ) {};

    bool build();

    ModelPtr getConnectionModel() const { return _model; };

    ModelPtr getParallelModel() const { return _pModel; };

    ParallelContactFEDescriptorPtr getParallelFiniteElementDescriptor() const { return _PFEDesc; };

    bool isParallel() const { return true; };
};

/**
 * @typedef ParallelContactNewPtr
 * @brief Pointeur intelligent vers un ParallelContactNew
 */
using ParallelContactNewPtr = std::shared_ptr< ParallelContactNew >;

class ParallelFrictionNew : public ParallelContactNew {

  public:
    /**
     * @typedef ParallelFrictionNewPtr
     * @brief Pointeur intelligent vers un ParallelFrictionNew
     */
    typedef std::shared_ptr< ParallelFrictionNew > ParallelFrictionNewPtr;

    /**
     * @brief Constructeur
     */
    ParallelFrictionNew() = delete;

    /**
     * @brief Constructeur
     */
    ParallelFrictionNew( const std::string name, const ModelPtr model, ParallelMeshPtr mesh )
        : ParallelContactNew( name, model, mesh, "CHAR_FROT" ) {};

    /**
     * @brief Constructeur
     */
    ParallelFrictionNew( const ModelPtr model, ParallelMeshPtr mesh )
        : ParallelFrictionNew( ResultNaming::getNewResultName(), model, mesh ) {};

    bool build() {
        AS_ASSERT( hasFriction() );

        return ParallelContactNew::build();
    };
};

/**
 * @typedef ParallelFrictionNewPtr
 * @brief Pointeur intelligent vers un ParallelFrictionNew
 */
using ParallelFrictionNewPtr = std::shared_ptr< ParallelFrictionNew >;

#endif /* ASTER_HAVE_MPI */
