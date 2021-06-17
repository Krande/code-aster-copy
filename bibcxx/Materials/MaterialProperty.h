#ifndef MATERIALBEHAVIOUR_H_
#define MATERIALBEHAVIOUR_H_

/**
 * @file MaterialProperty.h
 * @brief Fichier entete de la classe MaterialProperty
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

#include <iomanip>
#include <map>
#include <sstream>
#include <string>

#include "astercxx.h"
#include "aster_utils.h"
#include "Materials/BaseMaterialProperty.h"
#include "DataFields/Table.h"
#include "Functions/Formula.h"
#include "Functions/Function.h"
#include "Functions/Function2D.h"
#include "MemoryManager/JeveuxVector.h"


typedef std::vector< FunctionPtr > VectorFunction;


/**
 * @class MaterialProperty
 * @brief Classe fille de GenericMaterialProperty
 * @author Jean-Pierre Lefebvre
 */
class MaterialProperty : public GenericMaterialProperty {
    std::string capitalizeName( const std::string &nameInit ) {
        std::string name( nameInit );
        if ( !name.empty() ) {
            name[0] = std::toupper( name[0] );

            for ( std::size_t i = 1; i < name.length(); ++i )
                name[i] = std::tolower( name[i] );
        }
        return name;
    };

  public:
    /**
     * @brief Constructeur
     */
    MaterialProperty( const std::string asterName, const std::string asterNewName = "" )
        : GenericMaterialProperty( asterName, asterNewName ){};

    bool addPropertyReal( std::string name, const bool mandatory ) {
        return GenericMaterialProperty::addPropertyReal( capitalizeName( name ),
                                  ElementaryMaterialPropertyReal( name, mandatory ) );
    };

    bool addPropertyReal( std::string name, const ASTERDOUBLE &value, const bool mandatory ) {
        return GenericMaterialProperty::addPropertyReal( capitalizeName( name ),
                                  ElementaryMaterialPropertyReal( name, value, mandatory ) );
    };

    bool addPropertyComplex( std::string name, const bool mandatory ) {
        return GenericMaterialProperty::addPropertyComplex( capitalizeName( name ),
                                   ElementaryMaterialPropertyComplex( name, mandatory ) );
    };

    bool addPropertyString( std::string name, const bool mandatory ) {
        return GenericMaterialProperty::addPropertyString( capitalizeName( name ),
                                  ElementaryMaterialPropertyString( name, mandatory ) );
    };

    bool addPropertyString( std::string name, const std::string &value, const bool mandatory ) {
        return GenericMaterialProperty::addPropertyString( capitalizeName( name ),
                                  ElementaryMaterialPropertyString( name, value, mandatory ) );
    };

    bool addPropertyFunction( std::string name, const bool mandatory ) {
        return GenericMaterialProperty::addPropertyFunction( capitalizeName( name ),
                                    ElementaryMaterialPropertyDataStructure( name, mandatory ) );
    };

    bool addPropertyTable( std::string name, const bool mandatory ) {
        return GenericMaterialProperty::addPropertyTable( capitalizeName( name ),
                                 ElementaryMaterialPropertyTable( name, mandatory ) );
    };

    bool addPropertyVectorOfReal( std::string name, const bool mandatory ) {
        return GenericMaterialProperty::addPropertyVectorOfReal(
            capitalizeName( name ), ElementaryMaterialPropertyVectorReal( name, mandatory ) );
    };

    bool addPropertyVectorOfFunction( std::string name, const bool mandatory ) {
        return GenericMaterialProperty::addPropertyVectorOfFunction(
            capitalizeName( name ), ElementaryMaterialPropertyVectorFunction( name, mandatory ) );
    };

    /**
     * @brief Build ".RDEP"
     * @return true
     */
    bool computeTractionFunction( FunctionPtr &doubleValues ) const;

    /**
     * @brief Get name link to the class
     * @return name
     */
    std::string getName() { return _asterName; };
};

/** @typedef Pointeur intelligent vers un comportement materiau */
typedef boost::shared_ptr< MaterialProperty > MaterialPropertyPtr;


#endif
