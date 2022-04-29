#ifndef DISCRETEPROBLEM_H_
#define DISCRETEPROBLEM_H_

/**
 * @file DiscreteComputation.h
 * @brief Header of class DiscreteComputation
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

#include "astercxx.h"

#include "Behaviours/BehaviourProperty.h"
#include "Discretization/Calcul.h"
#include "LinearAlgebra/ElementaryMatrix.h"
#include "LinearAlgebra/ElementaryVector.h"
#include "Numbering/DOFNumbering.h"
#include "Studies/PhysicalProblem.h"

/**
 * @class DiscreteComputation
 * @brief Compute discrete operators (vectors and matrices)
 */
class DiscreteComputation {
  private:
    /** @brief Physical problem */
    PhysicalProblemPtr _phys_problem;

    /** @brief Compute B elementary matrices fo dualized boundary conditions */
    void baseDualStiffnessMatrix( CalculPtr &calcul,
                                  ElementaryMatrixDisplacementRealPtr &elemMatr );

    /** @brief Preparation for non-linear computations */
    CalculPtr createCalculForNonLinear( const std::string option,
                                        const ConstantFieldOnCellsRealPtr _timeFieldPrev,
                                        const ConstantFieldOnCellsRealPtr _timeFieldCurr,
                                        const FieldOnCellsRealPtr _externVarFieldPrev,
                                        const FieldOnCellsRealPtr _externVarFieldCurr,
                                        const VectorString &groupOfCells = VectorString() );

  public:
    /** @typedef DiscreteComputationPtr */
    typedef std::shared_ptr< DiscreteComputation > DiscreteComputationPtr;

    /** @brief Default constructor disabled */
    DiscreteComputation( void ) = delete;

    /**
     * @brief Constructor
     * @param PhysicalProblemPtr study
     */
    DiscreteComputation( const PhysicalProblemPtr &currPhysProblem )
        : _phys_problem( currPhysProblem ){};

    /** @brief Destructor */
    ~DiscreteComputation(){};

    /**
     * @brief Compute imposed displacement U_impo with Lagrange
     * @param time Time
     * @return Nodal field for imposed displacement
     */
    FieldOnNodesRealPtr imposedDisplacement( ASTERDOUBLE currTime = 0. );

    /**
     * @brief Compute Dirichlet reaction vector B^T * \lambda
     * @param time Time
     * @return Nodal field for Dirichlet reaction vector
     */
    FieldOnNodesRealPtr dualReaction( FieldOnNodesRealPtr lagr_curr );

    /**
     * @brief Compute Dirichlet imposed dualized displacement B * U
     * @param time Time
     * @return Nodal field for dualized displacement
     */
    FieldOnNodesRealPtr dualDisplacement( FieldOnNodesRealPtr disp_curr,
                                          ASTERDOUBLE scaling = 1.0 );

    /**
     * @brief Compute Neumann loads
     * @param TimeParameters Parameters for time
     * @return Nodal field for Neumann loads
     */
    FieldOnNodesRealPtr neumann( const VectorReal timeParameters,
                                 const FieldOnCellsRealPtr _externVarField = nullptr );

    /**
     * @brief Compute elementary matrices for mechanical stiffness (RIGI_MECA)
     * @param time Time
     * @return Elementary matrices for mechanical stiffness (RIGI_MECA)
     */
    ElementaryMatrixDisplacementRealPtr
    elasticStiffnessMatrix( const ASTERDOUBLE &time = 0.0, const ASTERINTEGER &modeFourier = 0,
                            const VectorString &groupOfCells = VectorString(),
                            const FieldOnCellsRealPtr _externVarField = nullptr );

    /**
     * @brief Compute elementary matrices for mass matrix (MASS_MECA)
     * @param time Time
     * @return Elementary matrices for mass matrix (MASS_MECA)
     */
    ElementaryMatrixDisplacementRealPtr massMatrix( ASTERDOUBLE time = 0. );

    /**
     * @brief Compute nodal field for kinematic boundary condition
     * @param time Time
     * @return Nodal field for kinematic boundary condition
     */
    FieldOnNodesRealPtr dirichletBC( const ASTERDOUBLE &time ) const;

    /**
     * @brief Compute nodal field for incremental kinematic boundary condition
     * @param time Time
     * @param disp_curr Current displacement
     * @return Nodal field for incremental kinematic boundary condition
     */
    FieldOnNodesRealPtr incrementalDirichletBC( const ASTERDOUBLE &time,
                                                const FieldOnNodesRealPtr disp_curr ) const;

    /**
     * @brief Compute B elementary matrices for dualized boundary conditions
     * @return Elementary matrices for dualized boundary conditions
     */
    ElementaryMatrixDisplacementRealPtr dualStiffnessMatrix();

    /**
     * @brief Get physical problem
     * @return Physical problem
     */
    PhysicalProblemPtr getPhysicalProblem() const { return _phys_problem; };

    /** @brief Create field for external state variables */
    FieldOnCellsRealPtr createExternalStateVariablesField( const ASTERDOUBLE time );

    /** @brief Create field for time */
    ConstantFieldOnCellsRealPtr createTimeField( const ASTERDOUBLE time );

    /** @brief Compute nodal field for external state variables RHS */
    FieldOnNodesRealPtr
    computeExternalStateVariablesLoad( const ASTERDOUBLE &time,
                                       const ConstantFieldOnCellsRealPtr _timeField,
                                       const FieldOnCellsRealPtr _externVarField ) const;

    /** @brief Compute field for external state variables reference values */
    FieldOnCellsRealPtr computeExternalStateVariablesReference() const;

    /**
     * @brief Compute internal forces, stress and internal state variables
     * @return Tuple with 5 objects:
     * field of exitcode
     * error code (integer)
     * internal state variables (VARI_ELGA)
     * Cauchy stress (SIEF_ELGA)
     * elementary vector of internal forces (`B^T \sigma`)
     */
    std::tuple< FieldOnCellsLongPtr, ASTERINTEGER, FieldOnCellsRealPtr, FieldOnCellsRealPtr,
                ElementaryVectorDisplacementRealPtr >
    computeInternalForces( const FieldOnNodesRealPtr displ, const FieldOnNodesRealPtr displ_incr,
                           const FieldOnCellsRealPtr stress, const FieldOnCellsRealPtr _internVar,
                           const ConstantFieldOnCellsRealPtr _timeFieldPrev,
                           const ConstantFieldOnCellsRealPtr _timeFieldCurr,
                           const VectorString &groupOfCells = VectorString() );

    /**
     * @brief Compute tangent matrix (not assembled)
     * @return Tuple with 5 objects:
     * field of exitcode
     * error code (integer)
     * internal state variables (VARI_ELGA)
     * Cauchy stress (SIEF_ELGA)
     * elementary vector of internal forces (`B^T \sigma`)
     * elementary tangent matrix
     */
    std::tuple< FieldOnCellsLongPtr, ASTERINTEGER, FieldOnCellsRealPtr, FieldOnCellsRealPtr,
                ElementaryVectorDisplacementRealPtr, ElementaryMatrixDisplacementRealPtr >
    computeTangentStiffnessMatrix( const FieldOnNodesRealPtr displ,
                                   const FieldOnNodesRealPtr displ_incr,
                                   const FieldOnCellsRealPtr stress,
                                   const FieldOnCellsRealPtr _internVar,
                                   const ConstantFieldOnCellsRealPtr _timeFieldPrev,
                                   const ConstantFieldOnCellsRealPtr _timeFieldCurr,
                                   const VectorString &groupOfCells = VectorString() );

    /**
     * @brief Compute tangent prediction matrix (not assembled)
     * @return Tuple with 4 objects:
     * field of maskcode
     * stress at prediction
     * elementary prediction matrix
     * elementary vector of internal forces (`B^T \sigma`)
     */
    std::tuple< FieldOnCellsLongPtr, FieldOnCellsRealPtr, ElementaryMatrixDisplacementRealPtr,
                ElementaryVectorDisplacementRealPtr >
    computeTangentPredictionMatrix( const FieldOnNodesRealPtr displ,
                                    const FieldOnNodesRealPtr displ_incr,
                                    const FieldOnCellsRealPtr stress,
                                    const FieldOnCellsRealPtr _internVar,
                                    const ConstantFieldOnCellsRealPtr _timeFieldPrev,
                                    const ConstantFieldOnCellsRealPtr _timeFieldCurr,
                                    const VectorString &groupOfCells = VectorString() );

    /**
     * @brief Compute elastic prediction matrix (not assembled)
     * @return Tuple with 1 objects:
     * elementary prediction matrix
     */
    std::tuple< ElementaryMatrixDisplacementRealPtr > computeElasticPredictionMatrix(
        const FieldOnNodesRealPtr displ, const FieldOnNodesRealPtr displ_incr,
        const FieldOnCellsRealPtr stress, const FieldOnCellsRealPtr _internVar,
        const ConstantFieldOnCellsRealPtr _timeFieldPrev,
        const ConstantFieldOnCellsRealPtr _timeFieldCurr );

    /**
     * @brief Compute secant prediction matrix (not assembled)
     * @return Tuple with 1 objects:
     * elementary prediction matrix
     */
    std::tuple< ElementaryMatrixDisplacementRealPtr > computeSecantPredictionMatrix(
        const FieldOnNodesRealPtr displ, const FieldOnNodesRealPtr displ_incr,
        const FieldOnCellsRealPtr stress, const FieldOnCellsRealPtr _internVar,
        const ConstantFieldOnCellsRealPtr _timeFieldPrev,
        const ConstantFieldOnCellsRealPtr _timeFieldCurr );
};

using DiscreteComputationPtr = std::shared_ptr< DiscreteComputation >;

#endif /* DISCRETEPROBLEM_H_ */
