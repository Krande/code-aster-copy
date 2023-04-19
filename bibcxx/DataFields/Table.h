#ifndef TABLE_H_
#define TABLE_H_

/**
 * @file Table.h
 * @brief Fichier entete de la classe Table
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
#ifdef ASTER_DEBUG_CXX
    bool _build_called = false;
#endif
    VectorString _parameters;
    std::map< std::string, JeveuxTypes > _typeByString;
    std::map< std::string, JeveuxVectorLong > _columnExists;
    std::map< std::string, JeveuxVectorLong > _columnLong;
    std::map< std::string, JeveuxVectorReal > _columnReal;
    std::map< std::string, JeveuxVectorComplex > _columnComplex;
    std::map< std::string, JeveuxVectorChar32 > _columnChar32;
    std::map< std::string, JeveuxVectorChar80 > _columnChar80;

  protected:
    std::map< std::string, JeveuxVectorChar8 > _columnChar8;
    std::map< std::string, JeveuxVectorChar16 > _columnChar16;
    std::map< std::string, JeveuxVectorChar24 > _columnChar24;
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
    Table( const std::string &, const std::string type = "TABLE" );

    /**
     * @brief Constructeur
     */
    Table();

    bool build();

    /**
     * @brief Return number of lines
     */
    int getNumberOfLines() const;

    /**
     * @brief Return parameters
     */
    const VectorString getParameters() const;

    /**
     * @brief Return Column type
     * @param column parameter
     */
    const std::string getColumnType( const std::string & ) const;

    /**
     * @brief Return Column values
     * @param column parameter
     */
    const std::tuple< VectorLong, VectorLong, VectorReal, VectorComplex, VectorString >
    getValues( const std::string & ) const;

    ~Table();
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
    TableOfFunctions( const std::string & );

    /**
     * @brief Constructeur
     */
    TableOfFunctions();

    /**
     * @brief Add function in TableOfFunctions
     * @param func function to add
     */
    void addFunction( GenericFunctionPtr );

    /**
     * @brief Get a function from his position
     * @param pos position
     */
    GenericFunctionPtr getFunction( int ) const;

    /**
     * @brief Get the number of functions referenced
     */
    int getNumberOfFunctions() const;
};

#endif /* TABLE_H_ */
