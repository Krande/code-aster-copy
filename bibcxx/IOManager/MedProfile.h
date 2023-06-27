#ifndef MEDPROFILE_H_
#define MEDPROFILE_H_

/**
 * @file MedProfile.h
 * @brief Fichier entete de la classe MedProfile
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

#include <iostream>
#include <memory>
#include <string>
#include <vector>

// aslint: disable=C3009

/**
 * @class MedProfile
 * @brief Med profile interface
 * @author Nicolas Sellenet
 */
class MedProfile {
  private:
    /** @brief entity id vector */
    std::vector< long > _indexVector;
    /** @brief profile name */
    const std::string _name;
    /** @brief profile size */
    const int _size;

  public:
    /**
     * @typedef MedProfilePtr
     * @brief Pointeur intelligent vers un MedProfile
     */
    typedef std::shared_ptr< MedProfile > MedProfilePtr;

    /**
     * @brief Constructeur
     */
    MedProfile( const std::string &name, int size ) : _name( name ), _size( size ) {};

    /** @brief get profile */
    const std::vector< long > &get() const { return _indexVector; };

    /** @brief get name */
    const std::string &getName() const { return _name; };

    /** @brief set profile */
    void set( const std::vector< long > &values ) { _indexVector = values; };
};

/**
 * @typedef MedProfilePtr
 * @brief Pointeur intelligent vers un MedProfile
 */
typedef std::shared_ptr< MedProfile > MedProfilePtr;

#endif /* MEDPROFILE_H_ */
