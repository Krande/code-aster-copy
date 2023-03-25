#pragma once

/**
 * @file Material.h
 * @brief Implementation of material datastructure.
 * @section LICENCE
 * Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
 * This file is part of code_aster.
 *
 * code_aster is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * code_aster is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "astercxx.h"

#include "aster_fort_material.h"

#include "DataStructures/DataStructure.h"
#include "Functions/Function.h"
#include "MemoryManager/JeveuxVector.h"
#include "Supervis/ResultNaming.h"

class MaterialListReal : public DataStructure {
  private:
    JeveuxVectorReal _list;

  public:
    MaterialListReal( const std::string &name, const VectorReal vect )
        : DataStructure( name, 16, "MATER_PROP_LISTR" ),
          _list( JeveuxVectorReal( getName() + ".LISV_R8", vect ) ) {};

    MaterialListReal( const std::string &name, const MaterialListReal &toCopy );
};

using MaterialListRealPtr = std::shared_ptr< MaterialListReal >;

class MaterialListFunc : public DataStructure {
  private:
    JeveuxVectorChar8 _list;

  public:
    MaterialListFunc( const std::string &name, const VectorString &vect );

    MaterialListFunc( const std::string &name, const MaterialListFunc &toCopy );
};

using MaterialListFuncPtr = std::shared_ptr< MaterialListFunc >;

class MaterialProperties : public DataStructure {
  private:
    JeveuxVectorReal _valR;
    JeveuxVectorComplex _valC;
    JeveuxVectorChar16 _valK;
    JeveuxVectorChar16 _ordr;
    JeveuxVectorLong _kord;

  public:
    MaterialProperties( const std::string &name );

    MaterialProperties( const std::string &name, const MaterialProperties &toCopy );

    MaterialProperties( const std::string &name, const int nbParam, const VectorReal valR,
                        const VectorComplex valC, const VectorString valK, const VectorString ordr,
                        const VectorLong kord );

    ASTERDOUBLE getValueReal( int idx );

    ASTERCOMPLEX getValueComplex( int idx );

    std::string getValueString( int idx );

    int getNumberOfReal() const { return _valR->size(); };

    int getNumberOfComplex() const { return _valC->size(); };

    int getNumberOfObjects() const {
        return ( _valK->size() - getNumberOfReal() - getNumberOfComplex() ) / 2;
    };
};

using MaterialPropertiesPtr = std::shared_ptr< MaterialProperties >;

class Material : public DataStructure {
  public:
    using MaterialPtr = std::shared_ptr< Material >;

  protected:
    /** @brief List of factor keywords / material names '.MATERIAU.NOMRC' */
    JeveuxVectorChar32 _names;
    /** @brief Function named '.&&RDEP' for traction function */
    FunctionPtr _rdep;

    std::vector< MaterialPropertiesPtr > _prop;
    /* List of .LISV_R8 vectors */
    std::vector< MaterialListRealPtr > _lisvR;
    /* List of .LISV_FO vectors */
    std::vector< MaterialListFuncPtr > _lisvF;
    /* List of the names of ListReal/ListFunc to simplify the creation of CodedMaterial */
    VectorString _nameList;

  private:
    bool build();
    std::string _cptName( int idx );
    MaterialPropertiesPtr matByName( std::string matName );
    int propIndex( std::string propName );

  public:
    Material() : Material( ResultNaming::getNewResultName() ) {};

    Material( const std::string &name );

    Material( const Material &toCopy );

    Material( const Material &toCopy, const VectorString MatToDelete );

    /**
     * @brief Add properties for a factor keyword / behaviour from user's keywords
     */
    void _addProperties( const std::string name, int nbParam, VectorReal valR, VectorComplex valC,
                         VectorString valK, VectorString ordr, VectorLong kord );

    ASTERINTEGER size();

    /**
     * @brief Store the LISTE_COEF values in an internal JeveuxVector
     */
    std::string _storeListReal( VectorReal vect );

    std::string _storeListFunc( VectorString vect );

    void _setTractionFunction( const std::string name, const std::string keyword,
                               GenericFunctionPtr &trac );

    VectorString getMaterialNames();

    ASTERDOUBLE getValueReal( const std::string matName, const std::string propName );

    ASTERCOMPLEX getValueComplex( const std::string matName, const std::string propName );

    GenericFunctionPtr getFunction( const std::string matName, const std::string propName );

    /* for compatibility */
    std::string getListName( const int position );
};

/**
 * @typedef MaterialPtr
 * @brief Smart pointer to a Material
 */
using MaterialPtr = std::shared_ptr< Material >;
using listOfMaterials = std::vector< MaterialPtr >;
