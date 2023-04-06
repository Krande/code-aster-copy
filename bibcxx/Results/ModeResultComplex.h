#ifndef MECHANICALMODECOMPLEXCONTAINER_H_
#define MECHANICALMODECOMPLEXCONTAINER_H_

/**
 * @file ModeResultComplex.h
 * @brief Fichier entete de la classe ModeResultComplex
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
#include "Results/ModeResult.h"
#include "Supervis/ResultNaming.h"

/**
 * @class ModeResultComplex
 * @brief Cette classe correspond a un mode_meca_c
 * On a choisi de définir un ModeResultComplex comme
   un résultat disposant, en plus des membres usuels d'un résultat, de champs aux noeuds complexes.
 * @author Natacha Béreux
 */
class ModeResultComplex : public ModeResult {
  private:
    using VectorOfComplexFieldsNodes = std::vector< FieldOnNodesComplexPtr >;
    using mapStrVOCFN = std::map< std::string, VectorOfComplexFieldsNodes >;
    using mapStrVOCFNIterator = mapStrVOCFN::iterator;
    using mapStrVOCFNValue = mapStrVOCFN::value_type;

    /** @brief Liste des champs aux noeuds */
    mapStrVOCFN _dictOfVectorOfComplexFieldsNodes;
    /** */
    StructureInterfacePtr _structureInterface;
    /** @brief Damping matrix */
    AssemblyMatrixDisplacementRealPtr _dampingMatrix;
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

  public:
    /**
     * @brief Constructeur
     */
    ModeResultComplex() : ModeResultComplex( ResultNaming::getNewResultName() ) {};

    /**
     * @brief Constructeur
     */
    ModeResultComplex( const std::string &name )
        : ModeResult( name, "MODE_MECA_C" ),
          _structureInterface( StructureInterfacePtr() ),
          _dampingMatrix( nullptr ),
          _rigidityDispDMatrix( nullptr ),
          _rigidityDispCMatrix( nullptr ),
          _rigidityTempDMatrix( nullptr ),
          _rigidityPressDMatrix( nullptr ),
          _rigidityGDMatrix( nullptr ),
          _rigidityGCMatrix( nullptr ) {};

    /**
     * @brief Set the damping matrix
     * @param matr AssemblyMatrixDisplacementRealPtr
     */
    bool setDampingMatrix( const AssemblyMatrixDisplacementRealPtr &matr ) {
        _dampingMatrix = matr;
        return true;
    };

    /**
     * @brief Set the rigidity matrix
     * @param matr AssemblyMatrixDisplacementRealPtr
     */
    bool setStiffnessMatrix( const AssemblyMatrixDisplacementRealPtr &matr ) {
        _rigidityDispDMatrix = matr;
        _rigidityDispCMatrix = nullptr;
        _rigidityTempDMatrix = nullptr;
        _rigidityPressDMatrix = nullptr;
        _rigidityGDMatrix = nullptr;
        _rigidityGCMatrix = nullptr;
        return true;
    };

    /**
     * @brief Set the rigidity matrix
     * @param matr AssemblyMatrixDisplacementComplexPtr
     */
    bool setStiffnessMatrix( const AssemblyMatrixDisplacementComplexPtr &matr ) {
        _rigidityDispDMatrix = nullptr;
        _rigidityDispCMatrix = matr;
        _rigidityTempDMatrix = nullptr;
        _rigidityPressDMatrix = nullptr;
        _rigidityGDMatrix = nullptr;
        _rigidityGCMatrix = nullptr;
        return true;
    };

    /**
     * @brief Set the rigidity matrix
     * @param matr AssemblyMatrixTemperatureRealPtr
     */
    bool setStiffnessMatrix( const AssemblyMatrixTemperatureRealPtr &matr ) {
        _rigidityDispDMatrix = nullptr;
        _rigidityDispCMatrix = nullptr;
        _rigidityTempDMatrix = matr;
        _rigidityPressDMatrix = nullptr;
        _rigidityGDMatrix = nullptr;
        _rigidityGCMatrix = nullptr;
        return true;
    };

    /**
     * @brief Set the rigidity matrix
     * @param matr AssemblyMatrixPressureRealPtr
     */
    bool setStiffnessMatrix( const AssemblyMatrixPressureRealPtr &matr ) {
        _rigidityDispDMatrix = nullptr;
        _rigidityDispCMatrix = nullptr;
        _rigidityTempDMatrix = nullptr;
        _rigidityPressDMatrix = matr;
        _rigidityGDMatrix = nullptr;
        _rigidityGCMatrix = nullptr;
        return true;
    };

    /**
     * @brief Set the rigidity matrix
     * @param matr GeneralizedAssemblyMatrixRealPtr
     */
    bool setStiffnessMatrix( const GeneralizedAssemblyMatrixRealPtr &matr ) {
        _rigidityDispDMatrix = nullptr;
        _rigidityDispCMatrix = nullptr;
        _rigidityTempDMatrix = nullptr;
        _rigidityPressDMatrix = nullptr;
        _rigidityGDMatrix = matr;
        _rigidityGCMatrix = nullptr;
        return true;
    };

    /**
     * @brief Set the rigidity matrix
     * @param matr GeneralizedAssemblyMatrixComplexPtr
     */
    bool setStiffnessMatrix( const GeneralizedAssemblyMatrixComplexPtr &matr ) {
        _rigidityDispDMatrix = nullptr;
        _rigidityDispCMatrix = nullptr;
        _rigidityTempDMatrix = nullptr;
        _rigidityPressDMatrix = nullptr;
        _rigidityGDMatrix = nullptr;
        _rigidityGCMatrix = matr;
        return true;
    };

    /**
     * @brief set interf_dyna
     * @param structureInterface objet StructureInterfacePtr
     */
    bool setStructureInterface( StructureInterfacePtr &structureInterface ) {
        _structureInterface = structureInterface;
        return true;
    };

    bool build() {
        BaseDOFNumberingPtr numeDdl( nullptr );
        if ( _dampingMatrix != nullptr )
            numeDdl = _dampingMatrix->getDOFNumbering();
        if ( _rigidityDispDMatrix != nullptr )
            numeDdl = _rigidityDispDMatrix->getDOFNumbering();
        if ( _rigidityDispCMatrix != nullptr )
            numeDdl = _rigidityDispCMatrix->getDOFNumbering();
        if ( _rigidityTempDMatrix != nullptr )
            numeDdl = _rigidityTempDMatrix->getDOFNumbering();
        if ( _rigidityPressDMatrix != nullptr )
            numeDdl = _rigidityPressDMatrix->getDOFNumbering();

        if ( numeDdl != nullptr ) {
            _mesh = numeDdl->getMesh();
        }
        return Result::build();
    };
};

using ModeResultComplexPtr = std::shared_ptr< ModeResultComplex >;

#endif /* MECHANICALMODECOMPLEXCONTAINER_H_ */
