/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

/**
 * @section LICENCE
 *   Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
 *   This file is part of code_aster.
 *
 *   code_aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   code_aster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 */

#include "astercxx.h"

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Numbering/GlobalEquationNumbering.h"

#pragma once

/**
 * @class ParallelGlobalEquationNumbering
 * @brief Class definissant un NUME_EQUA
 */
class ParallelGlobalEquationNumbering : public GlobalEquationNumbering {
  protected:
    /** @brief Objet Jeveux '.NULG' */
    JeveuxVectorLong _localToGlobal;
    /** @brief Objet Jeveux '.PDDL' */
    JeveuxVectorLong _localToRank;

  public:
    /**
     * @typedef GlobalEquationNumberingPtr
     * @brief Pointeur intelligent vers un ParallelGlobalEquationNumbering
     */
    typedef std::shared_ptr< ParallelGlobalEquationNumbering > ParallelGlobalEquationNumberingPtr;

    ParallelGlobalEquationNumbering( const std::string &baseName );

    /**
     * @brief Returns the vector of local to global numbering
     */
    const JeveuxVectorLong getLocalToGlobal() const { return _localToGlobal; };

    /**
     * @brief Returns the vector of the rank owning the local dof number
     */
    const JeveuxVectorLong getLocalToRank() const { return _localToRank; };

    bool isParallel() const { return true; };
};

/**
 * @typedef ParallelGlobalEquationNumberingPtr
 * @brief Enveloppe d'un pointeur intelligent vers un ParallelGlobalEquationNumbering
 */
typedef std::shared_ptr< ParallelGlobalEquationNumbering > ParallelGlobalEquationNumberingPtr;
