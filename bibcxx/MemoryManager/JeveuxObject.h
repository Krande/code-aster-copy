#ifndef JEVEUXOBJECT_H_
#define JEVEUXOBJECT_H_

/**
 * @file JeveuxObject.h
 * @brief Fichier entete de la classe JeveuxObject
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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

#include "aster_fort_jeveux.h"
#include "shared_vars.h"

#include <string>

/**
 * @enum JeveuxMemory
 * @brief To give memory space of Jeveux object (enum)
 */
enum JeveuxMemory { Permanent, Temporary };
/**
 * @def JeveuxMemoryTypesNames
 * @brief To give memory space of Jeveux object (string)
 */
static const std::string JeveuxMemoryTypesNames[2] = { "G", "V" };

/**
 * @def JeveuxMemoryTypesNames
 * @brief To give memory space of Jeveux object (string)
 */
static const int JeveuxNameMaxLength = 24;

/**
 * @class JeveuxObjectClass
 * @brief Cette classe permet de definir un objet Jeveux
 * @author Nicolas Sellenet
 */
class JeveuxObjectClass {
  protected:
    /** @brief Nom de l'objet Jeveux */
    std::string _name;
    /** @brief Mémoire d'allocation - Toujours globale en c++ */
    JeveuxMemory _mem;

  public:
    /**
     * @brief Constructeur
     * @param name Nom jeveux du vecteur
     */
    JeveuxObjectClass( const std::string &nom ) : _name( nom ), _mem( Permanent ) {};

    /**
     * @brief Destructeur
     */
    ~JeveuxObjectClass() {
        // #ifdef ASTER_DEBUG_CXX
        //         std::cout << "DEBUG: JeveuxObject.destr: " << _name << std::endl;
        // #endif
        if ( _name != "" && get_sh_jeveux_status() == 1 ) {
            CALLO_JEDETR( _name );
        }
    };

    bool exists() const;

    /**
     * @brief Return the name
     */
    std::string getName() const { return _name; };
};

#endif /* JEVEUXOBJECT_H_ */
