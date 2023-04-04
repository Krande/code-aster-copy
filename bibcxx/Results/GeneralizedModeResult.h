#ifndef GENERALIZEDMODERESULT_H_
#define GENERALIZEDMODERESULT_H_

/**
 * @file GeneralizedModeResult.h
 * @brief Fichier entete de la classe GeneralizedModeResult
 * @author Nicolas Tardieu
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

#include "LinearAlgebra/GeneralizedAssemblyMatrix.h"
#include "Modal/StaticMacroElement.h"
#include "Numbering/GeneralizedDOFNumbering.h"
#include "Results/FullResult.h"
#include "Supervis/ResultNaming.h"

/**
 * @class GeneralizedModeResult
 * @brief Cette classe correspond à un mode_gene
 * @author Nicolas Sellenet
 */
class GeneralizedModeResult : public FullResult {
  private:
    /** @brief Damping matrix */
    GeneralizedAssemblyMatrixRealPtr _dampingMatrix;
    /** @brief Stiffness double matrix */
    GeneralizedAssemblyMatrixRealPtr _rigidityRealMatrix;
    /** @brief Stiffness complex matrix */
    GeneralizedAssemblyMatrixComplexPtr _rigidityComplexMatrix;
    /** @brief generalized DOFNumbering */
    GeneralizedDOFNumberingPtr _genDOFNum;
    /** @brief Objet PROJ_MESU */
    ProjMesuPtr _projM;

  public:
    /**
     * @brief Constructeur
     * @todo  Ajouter les objets Jeveux de la SD
     */
    GeneralizedModeResult( const std::string &name )
        : FullResult( name, "MODE_GENE" ),
          _projM( new ProjMesu( ljust( getName(), 8 ) + ".PROJM" ) ),
          _rigidityRealMatrix( nullptr ),
          _rigidityComplexMatrix( nullptr ),
          _genDOFNum( nullptr ) {};

    /**
     * @brief Constructeur
     * @todo  Ajouter les objets Jeveux de la SD
     */
    GeneralizedModeResult() : GeneralizedModeResult( ResultNaming::getNewResultName() ) {};

    /**
     * @brief Get GeneralizedDOFNumering
     */
    GeneralizedDOFNumberingPtr getGeneralizedDOFNumbering() const {
        if ( _genDOFNum != nullptr )
            return _genDOFNum;
        throw std::runtime_error( "GeneralizedDOFNumbering is empty" );
    };

    /**
     * @brief Set the damping matrix
     * @param matr GeneralizedAssemblyMatrixRealPtr
     */
    bool setDampingMatrix( const GeneralizedAssemblyMatrixRealPtr &matr ) {
        _dampingMatrix = matr;
        return true;
    };

    /**
     * @brief Get the damping matrix
     * @param matr GeneralizedAssemblyMatrixRealPtr
     */
    GeneralizedAssemblyMatrixRealPtr getDampingMatrix( void ) const { return _dampingMatrix; };

    /**
     * @brief Set GeneralizedDOFNumering
     */
    bool setGeneralizedDOFNumbering( const GeneralizedDOFNumberingPtr &dofNum ) {
        if ( dofNum != nullptr ) {
            _genDOFNum = dofNum;
            //_fieldBuilder.addEquationNumbering( _genDOFNum->getDescription() );
            return true;
        }
        return false;
    };

    /**
     * @brief Set the rigidity matrix
     * @param matr GeneralizedAssemblyMatrixRealPtr
     */
    bool setStiffnessMatrix( const GeneralizedAssemblyMatrixRealPtr &matr ) {
        _rigidityComplexMatrix = nullptr;
        _rigidityRealMatrix = matr;
        return true;
    };

    /**
     * @brief Set the rigidity matrix
     * @param matr GeneralizedAssemblyMatrixComplexPtr
     */
    bool setStiffnessMatrix( const GeneralizedAssemblyMatrixComplexPtr &matr ) {
        _rigidityRealMatrix = nullptr;
        _rigidityComplexMatrix = matr;
        return true;
    };

    /**
     * @brief Get the stiffness matrix
     * @param matr GeneralizedAssemblyMatrixRealPtr
     */
    GeneralizedAssemblyMatrixRealPtr getGeneralizedStiffnessMatrixReal( void ) const {
        return _rigidityRealMatrix;
    };

    /**
     * @brief Get the stiffness matrix
     * @param matr GeneralizedAssemblyMatrixComplexPtr
     */
    GeneralizedAssemblyMatrixComplexPtr getGeneralizedStiffnessMatrixComplex( void ) const {
        return _rigidityComplexMatrix;
    };

    /**
     * @brief Get a generalized vector from its name and index
     * @param name name of the vector
     * @param index index number
     * @return GeneralizedAssemblyVectorRealPtr
     */
    GeneralizedAssemblyVectorRealPtr getGeneralizedVectorReal( const std::string name,
                                                               const ASTERINTEGER index ) const {
        return _dictOfMapOfGeneralizedVectorReal.at( name ).at( index );
    };

    /**
     * @brief Get a generalized vector from its name and index
     * @param name name of the vector
     * @param index index number
     * @return GeneralizedAssemblyVectorComplexPtr
     */
    GeneralizedAssemblyVectorComplexPtr
    getGeneralizedVectorComplex( const std::string name, const ASTERINTEGER index ) const {
        return _dictOfMapOfGeneralizedVectorComplex.at( name ).at( index );
    };

    bool build() { return Result::build(); };

    ~GeneralizedModeResult() {
        _fieldBuilder.clear();
        _dictOfMapOfFieldOnNodesReal.clear();
    };
};

/**
 * @typedef GeneralizedModeResultPtr
 * @brief Pointeur intelligent vers un GeneralizedModeResult
 */
typedef std::shared_ptr< GeneralizedModeResult > GeneralizedModeResultPtr;

#endif /* GENERALIZEDMODERESULT_H_ */
