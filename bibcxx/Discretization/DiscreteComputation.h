#ifndef DISCRETEPROBLEM_H_
#define DISCRETEPROBLEM_H_

/**
 * @file DiscreteComputation.h
 * @brief Header of class DiscreteComputation
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

#include "astercxx.h"

#include "Behaviours/BehaviourProperty.h"
#include "Discretization/Calcul.h"
#include "LinearAlgebra/ElementaryMatrix.h"
#include "LinearAlgebra/ElementaryVector.h"
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

    /** @brief Create time field */
    ConstantFieldOnCellsRealPtr createTimeField( const ASTERDOUBLE time,
                                                 const ASTERDOUBLE time_step = 0.0,
                                                 const ASTERDOUBLE theta = 0.0 ) const;

    /** @brief Create damping fluid field */
    ConstantFieldOnCellsLongPtr createDampingFluidField( const ASTERINTEGER damping,
                                                         const ASTERINTEGER onde_flui ) const;

    /** @brief Create wave type fluid field (for IMPE_MECA and AMOR_ACOU option) */
    ConstantFieldOnCellsLongPtr createWaveTypeFluidField( const ASTERINTEGER onde_flui ) const;

    /** @brief Preparation for non-linear computations */
    CalculPtr createCalculForNonLinear( const std::string option, const ASTERDOUBLE &time_prev,
                                        const ASTERDOUBLE &time_step,
                                        const FieldOnCellsRealPtr _externVarFieldPrev,
                                        const FieldOnCellsRealPtr _externVarFieldCurr,
                                        const VectorString &groupOfCells = VectorString() ) const;

    /** @brief Compute B elementary matrices fo dualized boundary conditions */
    void baseDualElasticStiffnessMatrix( CalculPtr &calcul,
                                         ElementaryMatrixDisplacementRealPtr &elemMatr ) const;

    /** @brief Compute B elementary matrices fo dualized thermal boundary conditions */
    void baseDualLinearConductivityMatrix( CalculPtr &calcul,
                                           ElementaryMatrixTemperatureRealPtr &elemMatr ) const;

    /** @brief Compute B elementary matrices fo dualized acoustic boundary conditions */
    void baseDualAcousticMatrix( CalculPtr &calcul,
                                 ElementaryMatrixPressureComplexPtr &elemMatr ) const;

    /** @brief Compute echange contributions to thermal matrix */
    void baseExchangeThermalMatrix( CalculPtr &calcul, ElementaryMatrixTemperatureRealPtr &elemMatr,
                                    const ASTERDOUBLE &time ) const;

    bool addTherImposedTerms( ElementaryVectorRealPtr elemVect, const ASTERDOUBLE time,
                              const ASTERDOUBLE time_step, const ASTERDOUBLE theta ) const;

    bool addMecaImposedTerms( ElementaryVectorRealPtr elemVect, const ASTERDOUBLE time ) const;

    bool addTherNeumannTerms( ElementaryVectorRealPtr elemVect, const ASTERDOUBLE time,
                              const ASTERDOUBLE time_step, const ASTERDOUBLE theta,
                              const FieldOnNodesRealPtr _previousNodalField = nullptr ) const;

    bool addMecaNeumannTerms( ElementaryVectorRealPtr elemVect, const ASTERDOUBLE time,
                              const ASTERDOUBLE time_step, const ASTERDOUBLE theta ) const;

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
        : _phys_problem( currPhysProblem ) {};

    /** @brief Destructor */
    ~DiscreteComputation() {};

    /** @brief Compute nodal field for external state variables RHS */
    FieldOnNodesRealPtr
    getExternalStateVariablesForces( const ASTERDOUBLE time,
                                     const FieldOnCellsRealPtr externVar = nullptr,
                                     const FieldOnCellsLongPtr maskField = nullptr ) const;

    /**
     * @brief Compute imposed displacement U_impo with Lagrange
     * @param time Time
     * @return Nodal field for imposed displacement
     */
    FieldOnNodesRealPtr getImposedDualBC( const ASTERDOUBLE time ) const {
        return this->getImposedDualBC( time, 0.0, 0.0 );
    };

    FieldOnNodesRealPtr getImposedDualBC( const ASTERDOUBLE time, const ASTERDOUBLE time_step,
                                          const ASTERDOUBLE theta ) const;
    /**
     * @brief Compute Dirichlet reaction vector B^T * \lambda
     * @param time time
     * @return Nodal field for Dirichlet reaction vector
     */
    FieldOnNodesRealPtr getDualForces( FieldOnNodesRealPtr lagr_curr ) const;

    /**
     * @brief Compute Dirichlet imposed dualized displacement B * U
     * @param time Time
     * @return Nodal field for dualized displacement
     */
    FieldOnNodesRealPtr getDualDisplacement( FieldOnNodesRealPtr disp_curr,
                                             ASTERDOUBLE scaling = 1.0 ) const;

    /**
     * @brief Compute Neumann loads
     * @param TimeParameters Parameters for time
     * @return Nodal field for Neumann loads
     */
    FieldOnNodesRealPtr
    getNeumannForces( const ASTERDOUBLE time = 0.0, const ASTERDOUBLE time_step = 0.0,
                      const ASTERDOUBLE theta = 1.0,
                      const FieldOnNodesRealPtr _previousPrimalField = nullptr ) const;

    /**
     * @brief Compute elementary matrices for mechanical stiffness (RIGI_MECA)
     */
    ElementaryMatrixDisplacementRealPtr
    getElasticStiffnessMatrix( const ASTERDOUBLE &time = 0.0, const ASTERINTEGER &modeFourier = 0,
                               const VectorString &groupOfCells = VectorString(),
                               const bool &with_dual = true ) const;

    /**
     * @brief Compute elementary matrices for mechanical stiffness (RIGI_FLUI_STRU)
     */
    ElementaryMatrixDisplacementRealPtr
    getFluidStructureStiffnessMatrix( const ASTERDOUBLE &time = 0.0,
                                      const ASTERINTEGER &modeFourier = 0,
                                      const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for rotational mechanical stiffness (RIGI_ROTA)
     */
    ElementaryMatrixDisplacementRealPtr getGeometricStiffnessMatrix(
        const FieldOnCellsRealPtr sief_elga, const FieldOnCellsRealPtr strx_elga = nullptr,
        const FieldOnNodesRealPtr displ = nullptr, const ASTERINTEGER &modeFourier = 0,
        const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for rotational mechanical stiffness (RIGI_ROTA)
     */
    ElementaryMatrixDisplacementRealPtr
    getRotationalStiffnessMatrix( const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for rotational mechanical stiffness (RIGI_GYRO)
     */
    ElementaryMatrixDisplacementRealPtr
    getGyroscopicStiffnessMatrix( const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for rotational mechanical damping (MECA_GYRO)
     */
    ElementaryMatrixDisplacementRealPtr
    getGyroscopicDampingMatrix( const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for thermal model (RIGI_THER)
     */
    ElementaryMatrixTemperatureRealPtr
    getLinearConductivityMatrix( const ASTERDOUBLE time, const ASTERINTEGER &modeFourier = -1,
                                 const VectorString &groupOfCells = VectorString(),
                                 const bool &with_dual = true ) const;

    /**
     * @brief Compute elementary matrices for thermal model (RIGI_THER_TANG)
     */
    ElementaryMatrixTemperatureRealPtr getTangentConductivityMatrix(
        const FieldOnNodesRealPtr temp, const FieldOnNodesRealPtr temp_step,
        const FieldOnCellsRealPtr &externVarCurr = nullptr,
        const VectorString &groupOfCells = VectorString(), const bool &with_dual = true ) const;

    /**
     * @brief Compute elementary matrices for acoustic model (RIGI_ACOU)
     */
    ElementaryMatrixPressureComplexPtr
    getLinearMobilityMatrix( const VectorString &groupOfCells = VectorString(),
                             const bool &with_dual = true ) const;

    /**
     * @brief Compute elementary matrices for mass matrix (MASS_MECA)
     */
    ElementaryMatrixDisplacementRealPtr
    getMechanicalMassMatrix( const bool diagonal, const ASTERDOUBLE &time = 0.0,
                             const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for mechanical stiffness (MASS_FLUI_STRU)
     */
    ElementaryMatrixDisplacementRealPtr
    getFluidStructureMassMatrix( const ASTERDOUBLE &time = 0.0,
                                 const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for mass matrix (MASS_THER)
     */
    ElementaryMatrixTemperatureRealPtr
    getLinearCapacityMatrix( const ASTERDOUBLE time,
                             const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for mass matrix (MASS_THER_TANG)
     */
    ElementaryMatrixTemperatureRealPtr
    getNonLinearCapacityMatrix( const FieldOnNodesRealPtr temp, const FieldOnNodesRealPtr temp_step,
                                const FieldOnCellsRealPtr &externVarCurr = nullptr,
                                const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for mass matrix (MASS_ACOU)
     */
    ElementaryMatrixPressureComplexPtr
    getCompressibilityMatrix( const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for damping matrix (AMOR_MECA)
     */
    ElementaryMatrixDisplacementRealPtr getMechanicalDampingMatrix(
        const ElementaryMatrixDisplacementRealPtr &massMatrix = nullptr,
        const ElementaryMatrixDisplacementRealPtr &stiffnessMatrix = nullptr,
        const ASTERDOUBLE &time = 0.0, const VectorString &groupOfCells = VectorString(),
        const ASTERINTEGER &flui_int = 1, const ASTERINTEGER &onde_flui = 1 ) const;
    /**
     * @brief Compute elementary matrices for damping matrix (AMOR_ACOU)
     */
    ElementaryMatrixPressureComplexPtr
    getImpedanceMatrix( const ASTERINTEGER &onde_flui = 1 ) const;

    /**
     * @brief Compute third order elementary matrices for absorbing fluid elements (IMPE_MECA)
     */
    ElementaryMatrixDisplacementRealPtr
    getImpedanceBoundaryMatrix( const VectorString &groupOfCells = VectorString(),
                                const ASTERINTEGER &onde_flui = 1 ) const;

    ElementaryMatrixDisplacementRealPtr
    getImpedanceWaveMatrix( const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for complex rigidity matrix (RIGI_MECA_HYST)
     */
    ElementaryMatrixDisplacementComplexPtr
    getHystereticStiffnessMatrix( const ElementaryMatrixDisplacementRealPtr &stiffnessMatrix,
                                  const ASTERDOUBLE &time = 0.0,
                                  const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute nodal field for kinematic boundary condition
     * @param time Time
     * @return Nodal field for kinematic boundary condition
     */
    template < typename T >
    std::shared_ptr< FieldOnNodes< T > > _getDirichletBC( const ASTERDOUBLE time = 0.0 ) const;
    FieldOnNodesRealPtr getMechanicalDirichletBC( const ASTERDOUBLE time = 0.0 ) const;
    FieldOnNodesRealPtr getThermalDirichletBC( const ASTERDOUBLE time = 0.0 ) const;
    FieldOnNodesComplexPtr getAcousticDirichletBC( const ASTERDOUBLE time = 0.0 ) const;

    /**
     * @brief Compute nodal field for incremental kinematic boundary condition
     * @param time Time
     * @param disp_curr Current displacement
     * @return Nodal field for incremental kinematic boundary condition
     */
    FieldOnNodesRealPtr getIncrementalDirichletBC( const ASTERDOUBLE &time,
                                                   const FieldOnNodesRealPtr disp_curr ) const;

    /**
     * @brief Compute B elementary matrices for dualized boundary conditions
     * @return Elementary matrices for dualized boundary conditions
     */
    ElementaryMatrixDisplacementRealPtr getDualElasticStiffnessMatrix() const;

    ElementaryMatrixPressureComplexPtr getDualLinearMobilityMatrix() const;

    ElementaryMatrixTemperatureRealPtr getDualLinearConductivityMatrix() const;

    ElementaryMatrixTemperatureRealPtr getExchangeThermalMatrix( const ASTERDOUBLE &time ) const;

    /**
     * @brief Get physical problem
     * @return Physical problem
     */
    PhysicalProblemPtr getPhysicalProblem() const { return _phys_problem; };

    FieldOnNodesRealPtr
    getTransientThermalForces( const ASTERDOUBLE time, const ASTERDOUBLE time_step,
                               const ASTERDOUBLE theta,
                               const FieldOnNodesRealPtr _previousPrimalField = nullptr ) const;

    FieldOnNodesRealPtr getNonLinearTransientThermalForces(
        const FieldOnNodesRealPtr temp, const FieldOnNodesRealPtr temp_step,
        const ASTERDOUBLE time_prev, const ASTERDOUBLE time_step, const ASTERDOUBLE theta,
        const FieldOnCellsRealPtr &externVarCurr = nullptr ) const;

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
    getInternalForces( const FieldOnNodesRealPtr displ, const FieldOnNodesRealPtr displ_step,
                       const FieldOnCellsRealPtr stress, const FieldOnCellsRealPtr internVar,
                       const ASTERDOUBLE &time_prev, const ASTERDOUBLE &time_step,
                       const FieldOnCellsRealPtr &externVarPrev = nullptr,
                       const FieldOnCellsRealPtr &externVarCurr = nullptr,
                       const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute internal forces, stress and internal state variables
     * @return Tuple with 5 objects:
     * field of internal forces (`B^T \sigma`)
     */
    FieldOnNodesRealPtr
    getInternalThermalForces( const FieldOnNodesRealPtr temp, const FieldOnNodesRealPtr temp_step,
                              const ASTERDOUBLE &time_prev, const ASTERDOUBLE &time_step,
                              const FieldOnCellsRealPtr &externVarCurr = nullptr,
                              const VectorString &groupOfCells = VectorString() ) const;

    // MASS_THER_RESI
    FieldOnNodesRealPtr
    getNonLinearCapacityForces( const FieldOnNodesRealPtr temp, const FieldOnNodesRealPtr temp_step,
                                const ASTERDOUBLE &time_prev, const ASTERDOUBLE &time_step,
                                const FieldOnCellsRealPtr &externVarCurr = nullptr,
                                const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute tangent matrix (not assembled)
     * @return Tuple with 3 objects:
     * field of exitcode
     * error code (integer)
     * elementary tangent matrix
     */
    std::tuple< FieldOnCellsLongPtr, ASTERINTEGER, ElementaryMatrixDisplacementRealPtr >
    getTangentStiffnessMatrix( const FieldOnNodesRealPtr displ,
                               const FieldOnNodesRealPtr displ_step,
                               const FieldOnCellsRealPtr stress,
                               const FieldOnCellsRealPtr internVar, const ASTERDOUBLE &time_prev,
                               const ASTERDOUBLE &time_step,
                               const FieldOnCellsRealPtr &externVarPrev = nullptr,
                               const FieldOnCellsRealPtr &externVarCurr = nullptr,
                               const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute tangent prediction matrix (not assembled)
     * @return Tuple with 3 objects:
     * field of exitcode
     * error code (integer)
     * elementary prediction matrix
     */
    std::tuple< FieldOnCellsLongPtr, ASTERINTEGER, ElementaryMatrixDisplacementRealPtr >
    getPredictionTangentStiffnessMatrix( const FieldOnNodesRealPtr displ,
                                         const FieldOnNodesRealPtr displ_step,
                                         const FieldOnCellsRealPtr stress,
                                         const FieldOnCellsRealPtr internVar,
                                         const ASTERDOUBLE &time_prev, const ASTERDOUBLE &time_step,
                                         const FieldOnCellsRealPtr &externVarPrev = nullptr,
                                         const FieldOnCellsRealPtr &externVarCurr = nullptr,
                                         const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute contact forces
     */
    FieldOnNodesRealPtr
    getContactForces( const MeshCoordinatesFieldPtr geom, const FieldOnNodesRealPtr displ,
                      const FieldOnNodesRealPtr displ_step, const ASTERDOUBLE &time_prev,
                      const ASTERDOUBLE &time_step, const FieldOnCellsRealPtr data,
                      const FieldOnNodesRealPtr coef_cont,
                      const FieldOnNodesRealPtr coef_frot ) const;

    /**
     * @brief Compute contact matrix
     */
    ElementaryMatrixDisplacementRealPtr
    getContactMatrix( const MeshCoordinatesFieldPtr geom, const FieldOnNodesRealPtr displ,
                      const FieldOnNodesRealPtr displ_step, const ASTERDOUBLE &time_prev,
                      const ASTERDOUBLE &time_step, const FieldOnCellsRealPtr data,
                      const FieldOnNodesRealPtr coef_cont,
                      const FieldOnNodesRealPtr coef_frot ) const;
};

using DiscreteComputationPtr = std::shared_ptr< DiscreteComputation >;

#endif /* DISCRETEPROBLEM_H_ */
