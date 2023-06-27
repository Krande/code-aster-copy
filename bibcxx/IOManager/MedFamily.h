#ifndef MEDFAMILY_H_
#define MEDFAMILY_H_

/**
 * @file MedFamily.h
 * @brief Fichier entete de la classe MedFamily
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

#include "astercxx.h"

#ifdef ASTER_HAVE_MED
#include "med.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

// aslint: disable=C3012

/**
 * @class MedFamily
 * @brief Med family interface
 * @author Nicolas Sellenet
 */
class MedFamily {
  private:
    const std::string _name;
    const med_int _num;
    const std::vector< std::string > _groups;

  public:
    /**
     * @typedef MedFamilyPtr
     * @brief Pointeur intelligent vers un MedFamily
     */
    typedef std::shared_ptr< MedFamily > MedFamilyPtr;

    /** @brief Constructor */
    MedFamily( const std::string &name, med_int num, const std::vector< std::string > &groups )
        : _name( name ), _num( num ), _groups( groups ) {};

    /** @brief Get group list for family */
    const std::vector< std::string > &getGroups() const { return _groups; };

    /** @brief Get family name */
    const std::string &getName() const { return _name; };

    /** @brief Get family id in med file */
    med_int getId() const { return _num; };
};

/**
 * @typedef MedFamilyPtr
 * @brief Pointeur intelligent vers un MedFamily
 */
typedef std::shared_ptr< MedFamily > MedFamilyPtr;

#endif
#endif /* MEDFAMILY_H_ */
