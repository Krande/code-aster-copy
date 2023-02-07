#ifndef TEMPORARYDATASTRUCTURENAMING_H_
#define TEMPORARYDATASTRUCTURENAMING_H_

/**
 * @file TemporaryDataStructureNaming.h
 * @brief Fichier entete de la classe TemporaryDataStructureNaming
 * @author Nicolas Sellenet
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
/* person_in_charge: nicolas.sellenet at edf.fr */

#include <sstream>
#include <string>

#include <assert.h>

/**
 * @class TemporaryDataStructureNaming
 * @brief Class to generate automatic JEVEUX name for a temporary object
 * @author Nicolas Sellenet
 */
class TemporaryDataStructureNaming {
  private:
    /** @brief Current index for automatic naming of Jeveux objects */
    static unsigned long int _number;

    /** @brief Maximum index for name of objects */
    static const unsigned long int maxNumberOfObjects = 268435455;

  public:
    /**
     * @brief Static member that returns a new name.
     */
    static std::string getNewTemporaryName( const int lengthName = 8 );
};

#endif /* TEMPORARYDATASTRUCTURENAMING_H_ */
