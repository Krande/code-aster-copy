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
 * @class MaterialPropertyClass
 * @brief Classe fille de GenericMaterialPropertyClass
 * @author Jean-Pierre Lefebvre
 */
class MaterialPropertyClass : public GenericMaterialPropertyClass {
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
    MaterialPropertyClass( const std::string asterName, const std::string asterNewName = "" )
        : GenericMaterialPropertyClass( asterName, asterNewName ){};

    bool addRealProperty( std::string name, const bool mandatory ) {
        return GenericMaterialPropertyClass::addRealProperty( capitalizeName( name ),
                                  ElementaryMaterialPropertyReal( name, mandatory ) );
    };

    bool addRealProperty( std::string name, const ASTERDOUBLE &value, const bool mandatory ) {
        return GenericMaterialPropertyClass::addRealProperty( capitalizeName( name ),
                                  ElementaryMaterialPropertyReal( name, value, mandatory ) );
    };

    bool addComplexProperty( std::string name, const bool mandatory ) {
        return GenericMaterialPropertyClass::addComplexProperty( capitalizeName( name ),
                                   ElementaryMaterialPropertyComplex( name, mandatory ) );
    };

    bool addStringProperty( std::string name, const bool mandatory ) {
        return GenericMaterialPropertyClass::addStringProperty( capitalizeName( name ),
                                  ElementaryMaterialPropertyString( name, mandatory ) );
    };

    bool addStringProperty( std::string name, const std::string &value, const bool mandatory ) {
        return GenericMaterialPropertyClass::addStringProperty( capitalizeName( name ),
                                  ElementaryMaterialPropertyString( name, value, mandatory ) );
    };

    bool addFunctionProperty( std::string name, const bool mandatory ) {
        return GenericMaterialPropertyClass::addFunctionProperty( capitalizeName( name ),
                                    ElementaryMaterialPropertyDataStructure( name, mandatory ) );
    };

    bool addTableProperty( std::string name, const bool mandatory ) {
        return GenericMaterialPropertyClass::addTableProperty( capitalizeName( name ),
                                 ElementaryMaterialPropertyTable( name, mandatory ) );
    };

    bool addVectorOfRealProperty( std::string name, const bool mandatory ) {
        return GenericMaterialPropertyClass::addVectorOfRealProperty(
            capitalizeName( name ), ElementaryMaterialPropertyVectorReal( name, mandatory ) );
    };

    bool addVectorOfFunctionProperty( std::string name, const bool mandatory ) {
        return GenericMaterialPropertyClass::addVectorOfFunctionProperty(
            capitalizeName( name ), ElementaryMaterialPropertyVectorFunction( name, mandatory ) );
    };

    /**
     * @brief Build ".RDEP"
     * @return true
     */
    bool buildTractionFunction( FunctionPtr &doubleValues ) const;

    /**
     * @brief Get name link to the class
     * @return name
     */
    std::string getName() { return _asterName; };
};

/** @typedef Pointeur intelligent vers un comportement materiau */
typedef boost::shared_ptr< MaterialPropertyClass > MaterialPropertyPtr;


#endif
