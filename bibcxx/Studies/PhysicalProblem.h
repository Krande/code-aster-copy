#ifndef PHYSICALPROBLEM_H_
#define PHYSICALPROBLEM_H_

/**
 * @file PhysicalProblem.h
 * @brief Fichier entete de la classe PhysicalProblem
 * @section LICENCE
 *   Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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

#include "aster_pybind.h"
#include "astercxx.h"

#include "Behaviours/BehaviourProperty.h"
#include "DataFields/FieldOnCells.h"
#include "Discretization/ElementaryCharacteristics.h"
#include "Loads/ListOfLoads.h"
#include "Materials/CodedMaterial.h"
#include "Materials/MaterialField.h"
//#include "Numbering/BaseDOFNumbering.h" not include because forward declaration

// Forward declaration
class BaseDOFNumbering;
using BaseDOFNumberingPtr = std::shared_ptr< BaseDOFNumbering >;

/**
 * @class PhysicalProblem
 * @brief Main class to describe physical problem
 */
class PhysicalProblem {
  private:
    /** @brief Mesh */
    BaseMeshPtr _mesh;

    /** @brief Model */
    ModelPtr _model;

    /** @brief Material field */
    MaterialFieldPtr _materialField;

    /** @brief List of loads */
    ListOfLoadsPtr _listOfLoads;

    /** @brief Elementary characteristics */
    ElementaryCharacteristicsPtr _elemChara;

    /** @brief Material parameters */
    CodedMaterialPtr _codedMater;

    /** @brief Behaviour properties */
    BehaviourPropertyPtr _behavProp;

    /** @brief Numbering */
    BaseDOFNumberingPtr _dofNume;

    /** @brief External state variable: reference field */
    FieldOnCellsRealPtr _externVarRefe;

  public:
    // No default constructor
    PhysicalProblem( void ) = delete;

    /** @brief Constructor */
    PhysicalProblem( const ModelPtr curModel, const MaterialFieldPtr curMat,
                     const ElementaryCharacteristicsPtr cara = nullptr );

    /** @brief Destructor */
    ~PhysicalProblem(){};

    /** @brief Add a load (mechanical or dirichlet) with function, formula */
    template < typename... Args >
    void addLoad( const Args &...a ) {
        _listOfLoads->addLoad( a... );
    };

    /** @brief Get material field */
    MaterialFieldPtr getMaterialField() const { return _materialField; };

    /** @brief Get model */
    ModelPtr getModel() const { return _model; };

    /** @brief Get mesh */
    BaseMeshPtr getMesh() const { return _mesh; };

    /** @brief Get elementary characteristics */
    ElementaryCharacteristicsPtr getElementaryCharacteristics() const { return _elemChara; };

    /** @brief Get material parameters*/
    CodedMaterialPtr getCodedMaterial() const { return _codedMater; };

    /** @brief Get list of loads */
    ListOfLoadsPtr getListOfLoads() const { return _listOfLoads; };

    /** @brief Get behaviour properties */
    BehaviourPropertyPtr getBehaviourProperty() const { return _behavProp; };

    /** @brief Get numbering of degrees of freedom */
    BaseDOFNumberingPtr getDOFNumbering() const { return _dofNume; };

    /** @brief Set numbering of degrees of freedom */
    void setDOFNumbering( const BaseDOFNumberingPtr dofNume );

    /** @brief Create numbering of degrees of freedom */
    bool computeDOFNumbering();

    /** @brief Create behaviour properties */
    void computeBehaviourProperty( py::object &keywords, const std::string &initialState = "NON",
                                   const ASTERINTEGER verbosity = 1 );

    /** @brief Create list of loads */
    bool computeListOfLoads() { return _listOfLoads->build( _model ); };

    /** @brief Compute field for external state variables reference values */
    void computeReferenceExternalStateVariables();

    /** @brief Get external state variables reference field */
    FieldOnCellsRealPtr getReferenceExternalStateVariables() const { return _externVarRefe; };

    /** @brief Get current external state variables reference field */
    FieldOnCellsRealPtr getExternalStateVariables( const ASTERDOUBLE &time ) const;
};

using PhysicalProblemPtr = std::shared_ptr< PhysicalProblem >;

#endif /* PHYSICALPROBLEM_H_ */
