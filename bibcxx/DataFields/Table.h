#ifndef TABLE_H_
#define TABLE_H_

/**
 * @file Table.h
 * @brief Fichier entete de la classe Table
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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

#include "DataStructures/DataStructure.h"
#include "Functions/Function.h"
#include "MemoryManager/JeveuxVector.h"
#include "Supervis/ResultNaming.h"

#include <string>

/**
 * @class Table
 * @brief Cette classe template permet de definir une table Aster
 * @author Nicolas Sellenet
 */
class Table : public DataStructure {
  protected:
    /** @brief Vecteur Jeveux '.TBBA' */
    JeveuxVectorChar8 _memoryLocation;
    /** @brief Vecteur Jeveux '.TBNP' */
    JeveuxVectorLong _description;
    /** @brief Vecteur Jeveux '.TBLP' */
    JeveuxVectorChar24 _parameterDescription;
    /** @brief Booléen indiquant si l'objet est vide */
    bool _isEmpty;

  public:
    /**
     * @typedef TablePtr
     * @brief Definition of a smart pointer to a Table
     */
    typedef std::shared_ptr< Table > TablePtr;

    // FIXME: Development documentation says 17 chars + "  ", for 'LG' logicals.

    /**
     * @brief Constructeur
     * @param name Nom Jeveux du champ aux noeuds
     */
    Table( const std::string &name, const std::string type = "TABLE" )
        : DataStructure( name, 19, type ),
          _memoryLocation( JeveuxVectorChar8( getName() + ".TBBA" ) ),
          _description( JeveuxVectorLong( getName() + ".TBNP" ) ),
          _parameterDescription( JeveuxVectorChar24( getName() + ".TBLP" ) ){};

    /**
     * @brief Constructeur
     */
    Table()
        : DataStructure( ResultNaming::getNewResultName(), 19, "TABLE" ),
          _memoryLocation( JeveuxVectorChar8( getName() + ".TBBA" ) ),
          _description( JeveuxVectorLong( getName() + ".TBNP" ) ),
          _parameterDescription( JeveuxVectorChar24( getName() + ".TBLP" ) ){};

    ~Table() {
        // #ifdef ASTER_DEBUG_CXX
        //         std::cout << "DEBUG: Table.destr: " << this->getName() << std::endl;
        // #endif
        if ( _parameterDescription->exists() && _description->exists() ) {
            _parameterDescription->updateValuePointer();
            _description->updateValuePointer();
            const int nbParam = ( *_description )[0];
            // Pour que la destruction soit effective, on crée des JeveuxVector pour déclencher
            // l'appel à JEDETR
            for ( int i = 0; i < nbParam; ++i ) {
                const JeveuxChar24 &name1 = ( *_parameterDescription )[i * 4 + 2];
                JeveuxVectorReal test1( name1.toString() );
                const JeveuxChar24 &name2 = ( *_parameterDescription )[i * 4 + 3];
                JeveuxVectorReal test2( name2.toString() );
            }
        }
    };
};

/**
 * @typedef TablePtr
 * @brief Definition of a smart pointer to a Table
 */
typedef std::shared_ptr< Table > TablePtr;

/**
 * @typedef TableOfFunctions
 * @brief Definition of TableOfFunctions (table_fonction)
 */
class TableOfFunctions : public Table {
  private:
    std::vector< GenericFunctionPtr > _vecOfFunctions;

  public:
    /**
     * @typedef TableOfFunctionsPtr
     * @brief Definition of a smart pointer to a TableOfFunctions
     */
    typedef std::shared_ptr< TableOfFunctions > TableOfFunctionsPtr;

    /**
     * @brief Constructeur
     * @param name Nom Jeveux du champ aux noeuds
     */
    TableOfFunctions( const std::string &name ) : Table( name, "TABLE_FONCTION" ){};

    /**
     * @brief Constructeur
     */
    TableOfFunctions() : Table( ResultNaming::getNewResultName(), "TABLE_FONCTION" ){};

    /**
     * @brief Add function in TableOfFunctions
     * @param func function to add
     */
    void addFunction( GenericFunctionPtr func ) { _vecOfFunctions.push_back( func ); };

    /**
     * @brief Get a function from his position
     * @param pos position
     */
    GenericFunctionPtr getFunction( int pos ) const {
        if ( pos < int( _vecOfFunctions.size() ) )
            return _vecOfFunctions[pos];
        return BaseFunctionPtr( nullptr );
    };

    /**
     * @brief Get the number of functions referenced
     */
    int getNumberOfFunctions() const { return _vecOfFunctions.size(); };
};

#endif /* TABLE_H_ */
