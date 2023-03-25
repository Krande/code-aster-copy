#ifndef MECHANICALMODERESULT_H_
#define MECHANICALMODERESULT_H_

/**
 * @file ModeResult.h
 * @brief Fichier entete de la classe ModeResult
 * @author Natacha Béreux
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

#include "Interfaces/StructureInterface.h"
#include "LinearAlgebra/AssemblyMatrix.h"
#include "LinearAlgebra/GeneralizedAssemblyMatrix.h"
#include "Numbering/DOFNumbering.h"
#include "Results/FullResult.h"
#include "Supervis/ResultNaming.h"

/**
 * @class ModeResult
 * @brief Cette classe correspond a un mode_meca
 * @author Nicolas Sellenet
 */
class ModeResult : public FullResult {
  private:
    StructureInterfacePtr _structureInterface;
    /** @brief Stiffness double displacement matrix */
    AssemblyMatrixDisplacementRealPtr _rigidityDispDMatrix;
    /** @brief Stiffness complex displacement matrix */
    AssemblyMatrixDisplacementComplexPtr _rigidityDispCMatrix;
    /** @brief Stiffness double temperature matrix */
    AssemblyMatrixTemperatureRealPtr _rigidityTempDMatrix;
    /** @brief Stiffness double pressure matrix */
    AssemblyMatrixPressureRealPtr _rigidityPressDMatrix;
    /** @brief Stiffness generalized double matrix */
    GeneralizedAssemblyMatrixRealPtr _rigidityGDMatrix;
    /** @brief Stiffness generalized complex matrix */
    GeneralizedAssemblyMatrixComplexPtr _rigidityGCMatrix;
    /** @brief Mass double displacement matrix */
    AssemblyMatrixDisplacementRealPtr _massDispDMatrix;
    /** @brief Mass complex displacement matrix */
    AssemblyMatrixDisplacementComplexPtr _massDispCMatrix;
    /** @brief Mass double temperature matrix */
    AssemblyMatrixTemperatureRealPtr _massTempDMatrix;
    /** @brief Mass double pressure matrix */
    AssemblyMatrixPressureRealPtr _massPressDMatrix;
    /** @brief Mass generalized double matrix */
    GeneralizedAssemblyMatrixRealPtr _massGDMatrix;
    /** @brief Mass generalized complex matrix */
    GeneralizedAssemblyMatrixComplexPtr _massGCMatrix;

  public:
    /**
     * @brief Constructeur
     */
    ModeResult() : ModeResult( ResultNaming::getNewResultName(), "MODE_MECA" ) {};

    /**
     * @brief Constructeur
     */
    ModeResult( const std::string &name, const std::string type = "MODE_MECA" )
        : FullResult( name, type ),
          _structureInterface( StructureInterfacePtr() ),
          _rigidityDispDMatrix( nullptr ),
          _rigidityDispCMatrix( nullptr ),
          _rigidityTempDMatrix( nullptr ),
          _rigidityPressDMatrix( nullptr ),
          _rigidityGDMatrix( nullptr ),
          _rigidityGCMatrix( nullptr ),
          _massDispDMatrix( nullptr ),
          _massDispCMatrix( nullptr ),
          _massTempDMatrix( nullptr ),
          _massPressDMatrix( nullptr ),
          _massGDMatrix( nullptr ),
          _massGCMatrix( nullptr ) {};

    /**
     * @brief Get the DOFNumbering
     */
    BaseDOFNumberingPtr getDOFNumbering() const {
        if ( _dofNum != nullptr )
            return _dofNum;
        if ( _rigidityDispDMatrix != nullptr )
            return _rigidityDispDMatrix->getDOFNumbering();
        if ( _rigidityTempDMatrix != nullptr )
            return _rigidityTempDMatrix->getDOFNumbering();
        return BaseDOFNumberingPtr( nullptr );
    };

    /**
     * @brief Get the rigidity matrix
     */
    AssemblyMatrixDisplacementComplexPtr getDisplacementComplexStiffnessMatrix() const {
        return _rigidityDispCMatrix;
    };

    /**
     * @brief Get the rigidity matrix
     */
    AssemblyMatrixDisplacementRealPtr getDisplacementRealStiffnessMatrix() const {
        return _rigidityDispDMatrix;
    };

    /**
     * @brief Get the rigidity matrix
     */
    AssemblyMatrixPressureRealPtr getPressureRealStiffnessMatrix() const {
        return _rigidityPressDMatrix;
    };

    /**
     * @brief Get the rigidity matrix
     */
    AssemblyMatrixTemperatureRealPtr getTemperatureRealStiffnessMatrix() const {
        return _rigidityTempDMatrix;
    };

    /**
     * @brief Set the rigidity matrix
     * @param matr AssemblyMatrixDisplacementRealPtr
     */
    void setStiffnessMatrix( const AssemblyMatrixDisplacementRealPtr &matr ) {
        _rigidityDispDMatrix = matr;
        _rigidityDispCMatrix = nullptr;
        _rigidityTempDMatrix = nullptr;
        _rigidityPressDMatrix = nullptr;
        _rigidityGDMatrix = nullptr;
        _rigidityGCMatrix = nullptr;
    };

    /**
     * @brief Set the rigidity matrix
     * @param matr AssemblyMatrixDisplacementComplexPtr
     */
    void setStiffnessMatrix( const AssemblyMatrixDisplacementComplexPtr &matr ) {
        _rigidityDispDMatrix = nullptr;
        _rigidityDispCMatrix = matr;
        _rigidityTempDMatrix = nullptr;
        _rigidityPressDMatrix = nullptr;
        _rigidityGDMatrix = nullptr;
        _rigidityGCMatrix = nullptr;
    };

    /**
     * @brief Set the rigidity matrix
     * @param matr AssemblyMatrixTemperatureRealPtr
     */
    void setStiffnessMatrix( const AssemblyMatrixTemperatureRealPtr &matr ) {
        _rigidityDispDMatrix = nullptr;
        _rigidityDispCMatrix = nullptr;
        _rigidityTempDMatrix = matr;
        _rigidityPressDMatrix = nullptr;
        _rigidityGDMatrix = nullptr;
        _rigidityGCMatrix = nullptr;
    };

    /**
     * @brief Set the rigidity matrix
     * @param matr AssemblyMatrixPressureRealPtr
     */
    void setStiffnessMatrix( const AssemblyMatrixPressureRealPtr &matr ) {
        _rigidityDispDMatrix = nullptr;
        _rigidityDispCMatrix = nullptr;
        _rigidityTempDMatrix = nullptr;
        _rigidityPressDMatrix = matr;
        _rigidityGDMatrix = nullptr;
        _rigidityGCMatrix = nullptr;
    };

    /**
     * @brief Set the rigidity matrix
     * @param matr GeneralizedAssemblyMatrixRealPtr
     */
    void setStiffnessMatrix( const GeneralizedAssemblyMatrixRealPtr &matr ) {
        _rigidityDispDMatrix = nullptr;
        _rigidityDispCMatrix = nullptr;
        _rigidityTempDMatrix = nullptr;
        _rigidityPressDMatrix = nullptr;
        _rigidityGDMatrix = matr;
        _rigidityGCMatrix = nullptr;
    };

    /**
     * @brief Set the rigidity matrix
     * @param matr GeneralizedAssemblyMatrixComplexPtr
     */
    void setStiffnessMatrix( const GeneralizedAssemblyMatrixComplexPtr &matr ) {
        _rigidityDispDMatrix = nullptr;
        _rigidityDispCMatrix = nullptr;
        _rigidityTempDMatrix = nullptr;
        _rigidityPressDMatrix = nullptr;
        _rigidityGDMatrix = nullptr;
        _rigidityGCMatrix = matr;
    };
    /**
     * @brief Get the mass matrix
     */
    AssemblyMatrixDisplacementComplexPtr getDisplacementComplexMassMatrix() const {
        return _massDispCMatrix;
    };

    /**
     * @brief Get the mass matrix
     */
    AssemblyMatrixDisplacementRealPtr getDisplacementRealMassMatrix() const {
        return _massDispDMatrix;
    };

    /**
     * @brief Get the mass matrix
     */
    AssemblyMatrixPressureRealPtr getPressureRealMassMatrix() const { return _massPressDMatrix; };

    /**
     * @brief Get the mass matrix
     */
    AssemblyMatrixTemperatureRealPtr getTemperatureRealMassMatrix() const {
        return _massTempDMatrix;
    };

    /**
     * @brief Set the mass matrix
     * @param matr AssemblyMatrixDisplacementRealPtr
     */
    void setMassMatrix( const AssemblyMatrixDisplacementRealPtr &matr ) {
        _massDispDMatrix = matr;
        _massDispCMatrix = nullptr;
        _massTempDMatrix = nullptr;
        _massPressDMatrix = nullptr;
        _massGDMatrix = nullptr;
        _massGCMatrix = nullptr;
    };

    /**
     * @brief Set the mass matrix
     * @param matr AssemblyMatrixDisplacementComplexPtr
     */
    void setMassMatrix( const AssemblyMatrixDisplacementComplexPtr &matr ) {
        _massDispDMatrix = nullptr;
        _massDispCMatrix = matr;
        _massTempDMatrix = nullptr;
        _massPressDMatrix = nullptr;
        _massGDMatrix = nullptr;
        _massGCMatrix = nullptr;
    };

    /**
     * @brief Set the mass matrix
     * @param matr AssemblyMatrixTemperatureRealPtr
     */
    void setMassMatrix( const AssemblyMatrixTemperatureRealPtr &matr ) {
        _massDispDMatrix = nullptr;
        _massDispCMatrix = nullptr;
        _massTempDMatrix = matr;
        _massPressDMatrix = nullptr;
        _massGDMatrix = nullptr;
        _massGCMatrix = nullptr;
    };

    /**
     * @brief Set the mass matrix
     * @param matr AssemblyMatrixPressureRealPtr
     */
    void setMassMatrix( const AssemblyMatrixPressureRealPtr &matr ) {
        _massDispDMatrix = nullptr;
        _massDispCMatrix = nullptr;
        _massTempDMatrix = nullptr;
        _massPressDMatrix = matr;
        _massGDMatrix = nullptr;
        _massGCMatrix = nullptr;
    };

    /**
     * @brief Set the mass matrix
     * @param matr GeneralizedAssemblyMatrixRealPtr
     */
    void setMassMatrix( const GeneralizedAssemblyMatrixRealPtr &matr ) {
        _massDispDMatrix = nullptr;
        _massDispCMatrix = nullptr;
        _massTempDMatrix = nullptr;
        _massPressDMatrix = nullptr;
        _massGDMatrix = matr;
        _massGCMatrix = nullptr;
    };

    /**
     * @brief Set the mass matrix
     * @param matr GeneralizedAssemblyMatrixComplexPtr
     */
    void setMassMatrix( const GeneralizedAssemblyMatrixComplexPtr &matr ) {
        _massDispDMatrix = nullptr;
        _massDispCMatrix = nullptr;
        _massTempDMatrix = nullptr;
        _massPressDMatrix = nullptr;
        _massGDMatrix = nullptr;
        _massGCMatrix = matr;
    };
    /**
     * @brief set interf_dyna
     * @param structureInterface objet StructureInterfacePtr
     */
    void setStructureInterface( StructureInterfacePtr &structureInterface ) {
        _structureInterface = structureInterface;
    };

    /**
     * @brief return number of dynamic modes
     */
    ASTERINTEGER getNumberOfDynamicModes() const {
        const std::string questi( "NB_MODES_DYN" );
        const std::string typeco( "RESULTAT" );
        ASTERINTEGER repi = 0, ier = 0;
        JeveuxChar32 repk( " " );
        const std::string arret( "F" );
        CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
        return repi;
    };

    /**
     * @brief return number of static modes
     */
    ASTERINTEGER getNumberOfStaticModes() const {
        const std::string questi( "NB_MODES_STA" );
        const std::string typeco( "RESULTAT" );
        ASTERINTEGER repi = 0, ier = 0;
        JeveuxChar32 repk( " " );
        const std::string arret( "F" );
        CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
        return repi;
    };

    bool build(
        const std::vector< FiniteElementDescriptorPtr > feds =
            std::vector< FiniteElementDescriptorPtr >(),
        const std::vector< EquationNumberingPtr > fnds = std::vector< EquationNumberingPtr >() ) {
        BaseDOFNumberingPtr numeDdl( nullptr );
        if ( _rigidityDispDMatrix != nullptr )
            numeDdl = _rigidityDispDMatrix->getDOFNumbering();
        if ( _rigidityDispCMatrix != nullptr )
            numeDdl = _rigidityDispCMatrix->getDOFNumbering();
        if ( _rigidityTempDMatrix != nullptr )
            numeDdl = _rigidityTempDMatrix->getDOFNumbering();
        if ( _rigidityPressDMatrix != nullptr )
            numeDdl = _rigidityPressDMatrix->getDOFNumbering();

        if ( numeDdl )
            setDOFNumbering( numeDdl );
        return Result::build( feds, fnds );
    };
};

/**
 * @typedef ModeResultPtr
 * @brief Pointeur intelligent vers un ModeResult
 */
typedef std::shared_ptr< ModeResult > ModeResultPtr;

#endif /* MECHANICALMODERESULT_H_ */
