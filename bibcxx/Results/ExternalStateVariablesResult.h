#ifndef EXTERNALVARIABLEEVOLUTIONCONTAINER_H_
#define EXTERNALVARIABLEEVOLUTIONCONTAINER_H_

/**
 * @file ExternalStateVariablesResult.h
 * @brief Fichier entete de la classe ExternalStateVariablesResult
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

/* person_in_charge: natacha.bereux at edf.fr */

#include "astercxx.h"

#include "Results/TransientResult.h"
#include "Supervis/ResultNaming.h"

/**
 * @class ExternalStateVariablesResult
 * @brief Cette classe correspond a un evol_varc, elle hérite de Result
          et stocke des champs
 * @author Natacha Béreux
 */
class ExternalStateVariablesResult : public TransientResult {
  private:
  public:
    /**
     * @brief Constructeur
     */
    ExternalStateVariablesResult( const std::string name = ResultNaming::getNewResultName() )
        : TransientResult( name, "EVOL_VARC" ) {};
};

/**
 * @typedef ExternalStateVariablesResultPtr
 * @brief Pointeur intelligent vers un ExternalStateVariablesResult
 */
typedef std::shared_ptr< ExternalStateVariablesResult > ExternalStateVariablesResultPtr;

#endif /* EXTERNALVARIABLEEVOLUTIONCONTAINER_H_ */
