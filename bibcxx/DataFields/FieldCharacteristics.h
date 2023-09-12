#ifndef FIELDCHARACTERISTICS_H_
#define FIELDCHARACTERISTICS_H_

/**
 * @file FieldCharacteristics.h
 * @brief Fichier entete de la classe FieldCharacteristics
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

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"

/**
 * @class FieldCharacteristics
 * @brief class which describe characteristics of a field from name (eg: DEPL)
 * @author Nicolas Sellenet
 */
class FieldCharacteristics {
  private:
    std::string _fieldQuantity;
    std::string _fieldSupport;
    std::string _option;
    std::string _parameter;
    std::string _name;

  public:
    /**
     * @typedef FieldCharacteristicsPtr
     * @brief Pointeur intelligent vers un FieldCharacteristics
     */
    typedef std::shared_ptr< FieldCharacteristics > FieldCharacteristicsPtr;

    /**
     * @brief Constructor
     * @param fieldName field name
     */
    FieldCharacteristics( const std::string &fieldName );

    const std::string &getQuantity() const { return _fieldQuantity; };

    const std::string &getLocalization() const { return _fieldSupport; };

    const std::string &getName() const { return _name; };

    const std::string &getOption() const { return _option; };

    const std::string &getParameter() const { return _parameter; };
};

#endif /* FIELDCHARACTERISTICS_H_ */
