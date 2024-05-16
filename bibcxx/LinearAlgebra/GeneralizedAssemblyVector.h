#ifndef GENERALIZEDASSEMBLYVECTOR_H_
#define GENERALIZEDASSEMBLYVECTOR_H_

/**
 * @file GeneralizedAssemblyVector.h
 * @brief Fichier entete de la classe GeneralizedAssemblyVector
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "Supervis/ResultNaming.h"

/**
 * @class GenericGeneralizedAssemblyVector
 * @brief Cette classe correspond a un vect_asse_gene
 * @author Nicolas Sellenet
 */
class GenericGeneralizedAssemblyVector : public DataStructure {
  private:
    /** @brief Objet Jeveux '.DESC' */
    JeveuxVectorLong _desc;
    /** @brief Objet Jeveux '.REFE' */
    JeveuxVectorChar24 _refe;

  public:
    /**
     * @brief Constructeur
     */
    GenericGeneralizedAssemblyVector( const std::string name )
        : DataStructure( name, 19, "VECT_ASSE_GENE" ),
          _desc( JeveuxVectorLong( getName() + ".DESC" ) ),
          _refe( JeveuxVectorChar24( getName() + ".REFE" ) ) {};
};

/**
 * @class GeneralizedAssemblyVector
 * @brief Cette classe correspond a un vect_asse_gene
 * @author Nicolas Sellenet
 */
template < class ValueType >
class GeneralizedAssemblyVector : public GenericGeneralizedAssemblyVector {
  private:
    /** @brief Objet Jeveux '.VALE' */
    JeveuxVector< ValueType > _vale;

  public:
    /**
     * @typedef GeneralizedAssemblyVectorPtr
     * @brief Pointeur intelligent vers un GeneralizedAssemblyVector
     */
    typedef std::shared_ptr< GeneralizedAssemblyVector< ValueType > > GeneralizedAssemblyVectorPtr;

    /**
     * @brief Constructeur
     */
    GeneralizedAssemblyVector() : GeneralizedAssemblyVector( ResultNaming::getNewResultName() ) {};

    /**
     * @brief Constructeur
     */

    GeneralizedAssemblyVector( const std::string name )
        : GenericGeneralizedAssemblyVector( name ),
          _vale( JeveuxVector< ValueType >( getName() + ".VALE" ) ) {};

    /**
     * @brief Get values of the field
     *
     */
    const JeveuxVector< ValueType > &getValues() const { return _vale; }

    void setValues( const std::vector< ValueType > &values ) {
        AS_ASSERT( values.size() == _vale->size() );
        *_vale = values;
    };
};

/** @typedef Definition d'une matrice assemblee généralisée de double */
typedef GeneralizedAssemblyVector< ASTERDOUBLE > GeneralizedAssemblyVectorReal;
/** @typedef Definition d'une matrice assemblee généralisée de complexe */
typedef GeneralizedAssemblyVector< ASTERCOMPLEX > GeneralizedAssemblyVectorComplex;

/**
 * @typedef GenericGeneralizedAssemblyVectorPtr
 * @brief Pointeur intelligent vers un GenericGeneralizedAssemblyVector
 */
typedef std::shared_ptr< GenericGeneralizedAssemblyVector > GenericGeneralizedAssemblyVectorPtr;

/**
 * @typedef GeneralizedAssemblyVectorRealPtr
 * @brief Pointeur intelligent vers un GeneralizedAssemblyVectorReal
 */
typedef std::shared_ptr< GeneralizedAssemblyVectorReal > GeneralizedAssemblyVectorRealPtr;

/**
 * @typedef GeneralizedAssemblyVectorComplexPtr
 * @brief Pointeur intelligent vers un GeneralizedAssemblyVectorComplex
 */
typedef std::shared_ptr< GeneralizedAssemblyVectorComplex > GeneralizedAssemblyVectorComplexPtr;

#endif /* GENERALIZEDASSEMBLYVECTOR_H_ */
