
#include "astercxx.h"

#ifdef ASTER_HAVE_MPI

#ifndef PARALLELDOFNUMBERING_H_
#define PARALLELDOFNUMBERING_H_

/**
 * @file ParallelDOFNumbering.h
 * @brief Fichier entete de la classe ParallelDOFNumbering
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "Numbering/BaseDOFNumbering.h"
#include "Numbering/ParallelEquationNumbering.h"

#include <unordered_map>

/**
 * @class ParallelDOFNumbering
 * @brief Class definissant un nume_ddl_p
 * @author Nicolas Sellenet
 */
class ParallelDOFNumbering : public BaseDOFNumbering {
  private:
    /** @brief Objet '.NUME' */
    ParallelEquationNumberingPtr _globalNumbering;

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

    ParallelDOFNumbering( const std::string name, const ParallelEquationNumberingPtr fdof,
                          const ModelPtr model );

    /**
     * @brief Constructeur
     * @param name nom souhaité de la sd (utile pour le BaseDOFNumbering d'une sd_resu)
     */
    ParallelDOFNumbering( const std::string &name );

    /**
     * @brief Returns the EquationNumberingPtr
     */
    EquationNumberingPtr getEquationNumbering() const { return _globalNumbering; };

    std::string getPhysicalQuantity() const { return _globalNumbering->getPhysicalQuantity(); };

    /**
     * @brief Methode permettant de savoir si l'objet est parallel
     * @return true
     */
    bool isParallel() { return true; };

    /**
     * @brief Are Lagrange Multipliers used for BC or MPC
     */
    bool useLagrangeDOF() const;

    /**
     * @brief Are Single Lagrange Multipliers used for BC or MPC
     */
    bool useSingleLagrangeDOF() const;

    /**
     * @brief Get The Component Associated To A Given Row
     */
    std::string getComponentFromDOF( const ASTERINTEGER dof, const bool local = false ) const;

    /**
     * @brief Get The Components Associated To A Given Node
     */
    VectorString getComponentFromNode( const ASTERINTEGER node, const bool local = false ) const;

    /**
     * @brief Get The Node Id Associated To A Given Row
     */
    ASTERINTEGER getNodeFromDOF( const ASTERINTEGER dof, const bool local = false ) const;

    /**
     * @brief Returns a vector with node index and component name for each DOFs
     */
    std::vector< std::pair< ASTERINTEGER, std::string > >
    getNodeAndComponentFromDOF( const bool local = false ) const;

    /**
     * @brief Returns a pair with node index and component name for given DOF
     */
    std::pair< ASTERINTEGER, std::string >
    getNodeAndComponentFromDOF( const ASTERINTEGER dof, const bool local = false ) const;

    /**
     * @brief Return true if a physical dof is Associated To A Given Row
     */
    bool isPhysicalDOF( const ASTERINTEGER dof, const bool local = false ) const;

    /**
     * @brief Get The total number of Dofs
     */
    ASTERINTEGER getNumberOfDOFs( const bool local = false ) const;

    /**
     * @brief Get Rows Associated to all Physical Dof
     */
    VectorLong getPhysicalDOFs( const bool local = false ) const;

    /**
     * @brief Get Rows Associated to all Ghost Dof
     */
    VectorLong getGhostDOFs( const bool local = false ) const;

    /**
     * @brief Get Rows owned locally (aka not Ghost)
     */
    VectorLong getNoGhostDOFs( const bool local = false ) const;

    /**
     * @brief Get Rows Associated to Lagrange Multipliers Dof
     */
    VectorLong getLagrangeDOFs( const bool local = false ) const;

    /**
     * @brief Get Rows Associated to Lagrange Multipliers Dof
     */
    std::map< ASTERINTEGER, VectorLong > getDictOfLagrangeDOFs( const bool local = false ) const;

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
     * @return Return the local number if the dof if present on the subdomain ; otherwise
     * raise an exception
     */
    ASTERINTEGER globalToLocalDOF( const ASTERINTEGER ) const;

    /**
     * @brief Return the global number of a local Dof
     * @return Return the global number if the dof if present on the subdomain ; otherwise
     * raise an exception
     */
    ASTERINTEGER localToGlobalDOF( const ASTERINTEGER ) const;

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

    void setMesh( const BaseMeshPtr mesh ) const { _globalNumbering->setMesh( mesh ); };

    /**
     * @brief Methode permettant de savoir si la numerotation est vide
     * @return true si la numerotation est vide
     */

    bool exists() const { return _globalNumbering->exists(); };

    bool build();
};

/**
 * @typedef ParallelDOFNumberingPtr
 * @brief Pointeur intelligent vers un ParallelDOFNumbering
 */
typedef std::shared_ptr< ParallelDOFNumbering > ParallelDOFNumberingPtr;

#endif /* PARALLELDOFNUMBERING_H_ */

#endif /* ASTER_HAVE_MPI */
