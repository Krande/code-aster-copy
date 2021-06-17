#ifndef MATERIALONMESHBUILDER_H_
#define MATERIALONMESHBUILDER_H_

/**
 * @file MaterialFieldBuilder.h
 * @brief Fichier entete de MaterialFieldBuilder
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "Materials/BaseExternalStateVariables.h"
#include "Materials/ListOfExternalStateVariables.h"
#include "Materials/MaterialField.h"
#include "astercxx.h"
#include <stdexcept>

/**
 * @class MaterialFieldBuilder
 * @author Nicolas Sellenet
 */
class MaterialFieldBuilder : public DataStructure {
    friend class MaterialField;

  protected:
    /**
     * @brief Build MaterialFieldPtr
     * @param curMater Material to build
     * @param curExternalStateVariable Input variables to add in MaterialFieldPtr
     */
    static void buildClass( MaterialField &curMater,
                        const ListOfExternalStateVariablesPtr &curExternalStateVariable = nullptr);

  public:
    /**
     * @typedef MaterialFieldBuilderPtr
     * @brief Pointeur intelligent vers un MaterialFieldBuilder
     */
    typedef boost::shared_ptr< MaterialFieldBuilder > MaterialFieldBuilderPtr;

    /**
     * @brief Build MaterialFieldPtr
     * @param curMater Material to build
     * @param curExternalStateVariable Input variables to add in MaterialFieldPtr
     */
    static MaterialFieldPtr build( MaterialFieldPtr &curMater,
                        const ListOfExternalStateVariablesPtr &curExternalStateVariable = nullptr);
};

#endif /* MATERIALONMESHBUILDER_H_ */
