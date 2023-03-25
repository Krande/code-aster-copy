#ifndef FULLRESULTSCONTAINER_H_
#define FULLRESULTSCONTAINER_H_

/**
 * @file FullResult.h
 * @brief Fichier entete de la classe FullResult
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

#include "Numbering/DOFNumbering.h"
#include "Results/DynamicResultsIndexing.h"
#include "Results/Result.h"
#include "Supervis/ResultNaming.h"

/**
 * @class FullResult
 * @brief Cette classe correspond à un sd_dyna_phys
 * @author Natacha Béreux
 */
class FullResult : public Result, DynamicResultsIndexing {
  protected:
    /** @brief the DOFNumbering */
    BaseDOFNumberingPtr _dofNum;

  public:
    /**
     * @brief Constructeur
     * @todo  Ajouter les objets Jeveux de la SD
     */
    FullResult( const std::string &name, const std::string &resuTyp )
        : Result( name, resuTyp ), DynamicResultsIndexing( getName() ), _dofNum( nullptr ) {};

    FullResult( const std::string &resuTyp )
        : FullResult( ResultNaming::getNewResultName(), resuTyp ) {};

    BaseDOFNumberingPtr getDOFNumbering() const { return _dofNum; };

    void setDOFNumbering( const BaseDOFNumberingPtr );
};

/**
 * @typedef FullResultPtr
 * @brief Pointeur intelligent vers un FullResult
 */
typedef std::shared_ptr< FullResult > FullResultPtr;

#endif /* FULLRESULTSCONTAINER_H_ */
