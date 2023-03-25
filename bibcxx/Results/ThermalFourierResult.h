#ifndef FOURIERTHERCONTAINER_H_
#define FOURIERTHERCONTAINER_H_

/**
 * @file ThermalFourierResult.h
 * @brief Fichier entete de la classe ThermalFourierResult
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

#include "Results/Result.h"
#include "Supervis/ResultNaming.h"

/**
 * @class ThermalFourierResult
 * @brief Cette classe correspond a un fourier_elas
 * @author Nicolas Sellenet
 */
class ThermalFourierResult : public Result {
  private:
  public:
    /**
     * @brief Constructeur
     */
    ThermalFourierResult( const std::string name = ResultNaming::getNewResultName() )
        : Result( name, "FOURIER_THER" ) {};
};

/**
 * @typedef ThermalFourierResultPtr
 * @brief Pointeur intelligent vers un ThermalFourierResult
 */
typedef std::shared_ptr< ThermalFourierResult > ThermalFourierResultPtr;

#endif /* FOURIERTHERCONTAINER_H_ */
