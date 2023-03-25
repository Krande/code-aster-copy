#ifndef ACOUSTICMODERESULT_H_
#define ACOUSTICMODERESULT_H_

/**
 * @file AcousticModeResult.h
 * @brief Fichier entete de la classe AcousticModeResult
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

#include "LinearAlgebra/AssemblyMatrix.h"
#include "Results/FullResult.h"
#include "Supervis/ResultNaming.h"

/**
 * @class AcousticModeResult
 * @brief Cette classe correspond a un mode_acou
 * @author Natacha Béreux
 */
class AcousticModeResult : public FullResult {
  private:
    /** @brief Stiffness displacement matrix */
    AssemblyMatrixPressureRealPtr _rigidityMatrix;

  public:
    /**
     * @brief Constructeur
     */
    AcousticModeResult( const std::string &name )
        : FullResult( name, "MODE_ACOU" ), _rigidityMatrix( nullptr ) {};

    AcousticModeResult() : AcousticModeResult( ResultNaming::getNewResultName() ) {};
    /**
     * @brief Set the rigidity matrix
     * @param matr AssemblyMatrixPressureRealPtr
     */
    bool setStiffnessMatrix( const AssemblyMatrixPressureRealPtr &matr ) {
        _rigidityMatrix = matr;
        return true;
    };

    bool build() { return Result::build(); };
};

/**
 * @typedef AcousticModeResultPtr
 * @brief Pointeur intelligent vers un AcousticModeResult
 */
typedef std::shared_ptr< AcousticModeResult > AcousticModeResultPtr;

#endif /* ACOUSTICMODERESULT_H_ */
