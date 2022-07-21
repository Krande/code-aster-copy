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
 * All methods are in the same header but implementation is splitted in several files
 * - DiscreteComputation.cxx for generic methods
 * - DiscreteComputationVector.cxx to compute vector thar are independent of physics
 * - DiscreteComputationMechanicalVector.cxx to compute vector for mechanics
 * - DiscreteComputationMechanicalMatrix.cxx to compute matrix for mechanics
 * - DiscreteComputationThermalMatrix.cxx to compute matrix for thermic
 * - DiscreteComputationThermalVector.cxx to compute vector for thermic
 */
class DiscreteComputation {
  private:
    /** @brief Physical problem */
    PhysicalProblemPtr _phys_problem;

    /** @brief Compute B elementary matrices fo dualized boundary conditions */
    void baseDualStiffnessMatrix( CalculPtr &calcul,
                                  ElementaryMatrixDisplacementRealPtr &elemMatr ) const;

    /** @brief Create time field */
    ConstantFieldOnCellsRealPtr createTimeField( const ASTERDOUBLE time_value,
                                                 const ASTERDOUBLE time_delta = 0.0,
                                                 const ASTERDOUBLE time_theta = 0.0 ) const;

    /** @brief Preparation for non-linear computations */
    CalculPtr createCalculForNonLinear( const std::string option, const ASTERDOUBLE &time_prev,
                                        const ASTERDOUBLE &time_step,
                                        const FieldOnCellsRealPtr _externVarFieldPrev,
                                        const FieldOnCellsRealPtr _externVarFieldCurr,
                                        const VectorString &groupOfCells = VectorString() ) const;

    /** @brief Compute B elementary matrices fo dualized thermal boundary conditions */
    void baseDualThermalMatrix( CalculPtr &calcul,
                                ElementaryMatrixTemperatureRealPtr &elemMatr ) const;

    /** @brief Compute echange contributions to thermal matrix */
    void baseExchangeThermalMatrix( CalculPtr &calcul, ElementaryMatrixTemperatureRealPtr &elemMatr,
                                    const ASTERDOUBLE time_value, const ASTERDOUBLE time_delta,
                                    const ASTERDOUBLE time_theta ) const;

    bool addTherImposedTerms( ElementaryVectorRealPtr elemVect, const ASTERDOUBLE time_value,
                              const ASTERDOUBLE time_delta, const ASTERDOUBLE time_theta ) const;

    bool addMecaImposedTerms( ElementaryVectorRealPtr elemVect,
                              const ASTERDOUBLE time_value ) const;

    bool addTherNeumannTerms( ElementaryVectorRealPtr elemVect, const ASTERDOUBLE time_value,
                              const ASTERDOUBLE time_delta, const ASTERDOUBLE time_theta,
                              const FieldOnCellsRealPtr _externVarField = nullptr,
                              const FieldOnNodesRealPtr _previousNodalField = nullptr ) const;

    bool addMecaNeumannTerms( ElementaryVectorRealPtr elemVect, const ASTERDOUBLE time_value,
                              const ASTERDOUBLE time_delta, const ASTERDOUBLE time_theta,
                              const FieldOnCellsRealPtr _externVarField = nullptr ) const;

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
    FieldOnNodesRealPtr imposedDualBC( const ASTERDOUBLE time_value ) const {
        return this->imposedDualBC( time_value, 0.0, 0.0 );
    };

    FieldOnNodesRealPtr imposedDualBC( const ASTERDOUBLE time_value, const ASTERDOUBLE time_delta,
                                       const ASTERDOUBLE time_theta ) const;
    /**
     * @brief Compute Dirichlet reaction vector B^T * \lambda
     * @param time Time_Value
     * @return Nodal field for Dirichlet reaction vector
     */
    FieldOnNodesRealPtr dualReaction( FieldOnNodesRealPtr lagr_curr ) const;

    /**
     * @brief Compute Dirichlet imposed dualized displacement B * U
     * @param time Time
     * @return Nodal field for dualized displacement
     */
    FieldOnNodesRealPtr dualDisplacement( FieldOnNodesRealPtr disp_curr,
                                          ASTERDOUBLE scaling = 1.0 ) const;

    /**
     * @brief Compute Neumann loads
     * @param TimeParameters Parameters for time
     * @return Nodal field for Neumann loads
     */
    FieldOnNodesRealPtr neumann( const ASTERDOUBLE time_value, const ASTERDOUBLE time_delta,
                                 const ASTERDOUBLE time_theta,
                                 const FieldOnCellsRealPtr _externVarField = nullptr,
                                 const FieldOnNodesRealPtr _previousPrimalField = nullptr ) const;
    /**
     * @brief Compute elementary matrices for mechanical stiffness (RIGI_MECA)
     * @param time_value Time
     * @param groupofCells GROUP_MA
     * @return Elementary matrices for mechanical stiffness (RIGI_MECA)
     */
    ElementaryMatrixDisplacementRealPtr
    elasticStiffnessMatrix( const ASTERDOUBLE &time_value = 0.0,
                            const ASTERINTEGER &modeFourier = 0,
                            const VectorString &groupOfCells = VectorString(),
                            const FieldOnCellsRealPtr _externVarField = nullptr ) const;

    /**
     * @brief Compute elementary matrices for thermal model (RIGI_THER)
     * @param time_value Time
     * @param groupofCells GROUP_MA
     * @return Elementary matrices for thermal model (RIGI_THER)
     */
    ElementaryMatrixTemperatureRealPtr
    linearConductivityMatrix( const ASTERDOUBLE time_value, const ASTERDOUBLE time_delta,
                              const ASTERDOUBLE time_theta, const ASTERINTEGER &modeFourier = -1,
                              const VectorString &groupOfCells = VectorString(),
                              const FieldOnCellsRealPtr _externVarField = nullptr ) const;

    /**
     * @brief Compute elementary matrices for mass matrix (MASS_MECA)
     * @param time Time
     * @param groupofCells GROUP_MA
     * @return Elementary matrices for mass matrix (MASS_MECA)
     */
    ElementaryMatrixDisplacementRealPtr
    massMatrix( const ASTERDOUBLE &time_value = 0.0,
                const VectorString &groupOfCells = VectorString(),
                const FieldOnCellsRealPtr _externVarField = nullptr ) const;

    /**
     * @brief Compute elementary matrices for mass matrix (MASS_THER)
     * @param time Time
     * @param groupofCells GROUP_MA
     * @return Elementary matrices for mass matrix (MASS_THER)
     */
    ElementaryMatrixTemperatureRealPtr
    linearCapacityMatrix( const ASTERDOUBLE time_value, const ASTERDOUBLE time_delta,
                          const ASTERDOUBLE time_theta,
                          const VectorString &groupOfCells = VectorString(),
                          const FieldOnCellsRealPtr _externVarField = nullptr ) const;

    /**
     * @brief Compute elementary matrices for damping matrix (AMOR_MECA)
     * @param massMatrix Elementary matrices for mass matrix
     * @param stiffnessMatrix  Elementary matrices for mechanical stiffness
     * @return Elementary matrices for damping matrix (AMOR_MECA)
     */
    ElementaryMatrixDisplacementRealPtr
    dampingMatrix( const ElementaryMatrixDisplacementRealPtr &massMatrix = nullptr,
                   const ElementaryMatrixDisplacementRealPtr &stiffnessMatrix = nullptr,
                   const ASTERDOUBLE &time_value = 0.0,
                   const VectorString &groupOfCells = VectorString(),
                   const FieldOnCellsRealPtr _externVarField = nullptr ) const;

    /**
     * @brief Compute nodal field for kinematic boundary condition
     * @param time Time
     * @return Nodal field for kinematic boundary condition
     */
    FieldOnNodesRealPtr dirichletBC( const ASTERDOUBLE &time_value ) const;

    /**
     * @brief Compute nodal field for incremental kinematic boundary condition
     * @param time Time
     * @param disp_curr Current displacement
     * @return Nodal field for incremental kinematic boundary condition
     */
    FieldOnNodesRealPtr incrementalDirichletBC( const ASTERDOUBLE &time_value,
                                                const FieldOnNodesRealPtr disp_curr ) const;

    /**
     * @brief Compute B elementary matrices for dualized boundary conditions
     * @return Elementary matrices for dualized boundary conditions
     */
    ElementaryMatrixDisplacementRealPtr dualStiffnessMatrix() const;

    /**
     * @brief Get physical problem
     * @return Physical problem
     */
    PhysicalProblemPtr getPhysicalProblem() const { return _phys_problem; };

    /** @brief Create field for external state variables */
    FieldOnCellsRealPtr createExternalStateVariablesField( const ASTERDOUBLE time_value ) const;

    /** @brief Compute nodal field for external state variables RHS */
    FieldOnNodesRealPtr
    computeExternalStateVariablesLoad( const ASTERDOUBLE time_value,
                                       const FieldOnCellsRealPtr _externVarField ) const;

    /** @brief Compute nodal field for external state variables RHS */
    FieldOnNodesRealPtr
    transientThermalLoad( const ASTERDOUBLE time_value, const ASTERDOUBLE time_delta,
                          const ASTERDOUBLE time_theta, const FieldOnCellsRealPtr _externVarField,
                          const FieldOnNodesRealPtr _previousPrimalField = nullptr ) const;

    /** @brief Compute field for external state variables reference values */
    FieldOnCellsRealPtr computeExternalStateVariablesReference() const;

    /**
     * @brief Compute internal forces, stress and internal state variables
     * @return Tuple with 5 objects:
     * field of exitcode
     * error code (integer)
     * internal state variables (VARI_ELGA)
     * Cauchy stress (SIEF_ELGA)
     * field of internal forces (`B^T \sigma`)
     */
    std::tuple< FieldOnCellsLongPtr, ASTERINTEGER, FieldOnCellsRealPtr, FieldOnCellsRealPtr,
                FieldOnNodesRealPtr >
    computeInternalForces( const FieldOnNodesRealPtr displ, const FieldOnNodesRealPtr displ_step,
                           const FieldOnCellsRealPtr stress, const FieldOnCellsRealPtr internVar,
                           const ASTERDOUBLE &time_prev, const ASTERDOUBLE &time_step,
                           const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute tangent matrix (not assembled)
     * @return Tuple with 3 objects:
     * field of exitcode
     * error code (integer)
     * elementary tangent matrix
     */
    std::tuple< FieldOnCellsLongPtr, ASTERINTEGER, ElementaryMatrixDisplacementRealPtr >
    computeTangentStiffnessMatrix( const FieldOnNodesRealPtr displ,
                                   const FieldOnNodesRealPtr displ_step,
                                   const FieldOnCellsRealPtr stress,
                                   const FieldOnCellsRealPtr internVar,
                                   const ASTERDOUBLE &time_prev, const ASTERDOUBLE &time_step,
                                   const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute tangent prediction matrix (not assembled)
     * @return Tuple with 3 objects:
     * field of exitcode
     * error code (integer)
     * elementary prediction matrix
     */
    std::tuple< FieldOnCellsLongPtr, ASTERINTEGER, ElementaryMatrixDisplacementRealPtr >
    computeTangentPredictionMatrix( const FieldOnNodesRealPtr displ,
                                    const FieldOnNodesRealPtr displ_step,
                                    const FieldOnCellsRealPtr stress,
                                    const FieldOnCellsRealPtr internVar,
                                    const ASTERDOUBLE &time_prev, const ASTERDOUBLE &time_step,
                                    const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute contact forces
     */
    FieldOnNodesRealPtr contactForces( const MeshCoordinatesFieldPtr geom,
                                       const FieldOnNodesRealPtr displ,
                                       const FieldOnNodesRealPtr displ_step,
                                       const FieldOnCellsRealPtr data ) const;

    /**
     * @brief Compute contact matrix
     */
    ElementaryMatrixDisplacementRealPtr contactMatrix( const MeshCoordinatesFieldPtr geom,
                                                       const FieldOnNodesRealPtr displ,
                                                       const FieldOnNodesRealPtr displ_step,
                                                       const FieldOnCellsRealPtr data ) const;
};

using DiscreteComputationPtr = std::shared_ptr< DiscreteComputation >;

#endif /* DISCRETEPROBLEM_H_ */
