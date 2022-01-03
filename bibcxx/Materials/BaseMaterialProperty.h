#ifndef BASEMATERIALBEHAVIOUR_H_
#define BASEMATERIALBEHAVIOUR_H_

/**
 * @file BaseMaterialProperty.h
 * @brief Fichier entete de la classe BaseMaterialProperty
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

#include <iomanip>
#include <map>
#include <sstream>
#include <string>

#include "DataFields/Table.h"
#include "Functions/Formula.h"
#include "Functions/Function.h"
#include "Functions/Function2D.h"
#include "MemoryManager/JeveuxVector.h"
#include "Utilities/ConvertibleValue.h"
#include "aster_utils.h"
#include "astercxx.h"

typedef std::vector< FunctionPtr > VectorFunction;



typedef ConvertibleValue< std::string, ASTERDOUBLE > StringToValueReal;

/**
 * @struct template AllowedMaterialPropertyType
 * @brief Structure permettant de limiter le type instanciable de MaterialProperty
 * @author Nicolas Sellenet
 */
template < typename T > struct AllowedMaterialPropertyType;

template <> struct AllowedMaterialPropertyType< ASTERDOUBLE > {};

template <> struct AllowedMaterialPropertyType< ASTERCOMPLEX > {};

template <> struct AllowedMaterialPropertyType< std::string > {};

template <> struct AllowedMaterialPropertyType< FunctionPtr > {};

template <> struct AllowedMaterialPropertyType< TablePtr > {};

template <> struct AllowedMaterialPropertyType< Function2DPtr > {};

template <> struct AllowedMaterialPropertyType< FormulaPtr > {};

template <> struct AllowedMaterialPropertyType< GenericFunctionPtr > {};

template <> struct AllowedMaterialPropertyType< VectorReal > {};

template <> struct AllowedMaterialPropertyType< VectorFunction > {};

template <> struct AllowedMaterialPropertyType< StringToValueReal > {};

class GenericMaterialProperty;

template < typename T1 > struct is_convertible;

template < class T > struct is_convertible {
    typedef T value_type;
    typedef T init_value;
};

template <> struct is_convertible< StringToValueReal > {
    typedef typename StringToValueReal::ReturnValue value_type;
    typedef typename StringToValueReal::BaseValue init_value;
};

/**
 * @class ElementaryMaterialProperty
 * @brief Cette classe template permet de definir un type elementaire de propriete materielle
 * @author Nicolas Sellenet
 * @todo on pourrait detemplatiser cette classe pour qu'elle prenne soit des doubles soit des fct
 *       on pourrait alors fusionner elas et elasFo par exemple
 */
template < class ValueType >
class ElementaryMaterialProperty : private AllowedMaterialPropertyType< ValueType > {
  public:
    typedef typename is_convertible< ValueType >::value_type ReturnValue;
    typedef typename is_convertible< ValueType >::init_value BaseValue;

  protected:
    /** @brief Nom Aster du type elementaire de propriete materielle */
    // ex : "NU" pour le coefficient de Poisson
    std::string _name;
    /** @brief Booleen qui precise si la propriété est obligatoire */
    bool _isMandatory;
    /** @brief Description de parametre, ex : "Young's modulus" */
    std::string _description;
    /** @brief Valeur du parametre (double, FunctionPtr, ...) */
    ValueType _value;
    /** @brief Booleen qui precise si la propriété a été initialisée */
    bool _existsValue;

  public:
    /**
     * @brief Constructeur vide (utile pour ajouter une ElementaryMaterialProperty a une
     std::map
     */
    ElementaryMaterialProperty(){};

    /**
     * @brief Constructeur
     * @param name Nom Aster du parametre materiau (ex : "NU")
     * @param description Description libre
     */
    ElementaryMaterialProperty( const std::string name, const bool isMandatory )
        : _name( name ), _isMandatory( isMandatory ), _description( "" ), _existsValue( false ){};

    /**
     * @brief Constructeur
     * @param name Nom Aster du parametre materiau (ex : "NU")
     * @param ValueType Valeur par défaut
     * @param description Description libre
     */
    ElementaryMaterialProperty( const std::string name, const ValueType &currentValue,
                              const bool isMandatory )
        : _name( name ), _isMandatory( isMandatory ), _description( "" ), _value( currentValue ),
          _existsValue( true ){};

    /**
     * @brief Recuperation de la valeur du parametre
     * @return le nom Aster du parametre
     */
    const std::string &getName() const { return _name; };

    /**
     * @brief Recuperation de la valeur du parametre
     * @return la valeur du parametre
     */
    template < typename T = ValueType >
    typename std::enable_if< std::is_same< T, StringToValueReal >::value,
                             const ReturnValue & >::type
    getValue() const {
        return _value.getValue();
    };

    /**
     * @brief Recuperation de la valeur du parametre
     * @return la valeur du parametre
     */
    template < typename T = ValueType >
    typename std::enable_if< !std::is_same< T, StringToValueReal >::value,
                             const ReturnValue & >::type
    getValue() const {
        return _value;
    };

    /**
     * @brief Cette propriété est-elle obligatoire ?
     * @return true si la propriété est obligatoire
     */
    bool isMandatory() const { return _isMandatory; };

    /**
     * @brief Cette propriété a-t-elle une valeur fixée par l'utilisateur ?
     * @return true si la valeur a été précisée
     */
    template < typename T = ValueType >
    typename std::enable_if< std::is_same< T, StringToValueReal >::value, bool >::type
    hasValue() const {
        return _value.hasValue();
    };

    /**
     * @brief Cette propriété a-t-elle une valeur fixée par l'utilisateur ?
     * @return true si la valeur a été précisée
     */
    template < typename T = ValueType >
    typename std::enable_if< !std::is_same< T, StringToValueReal >::value, bool >::type
    hasValue() const {
        return _existsValue;
    };

    /**
     * @brief Fonction servant a fixer la valeur du parametre
     * @param currentValue valeur donnee par l'utilisateur
     */
    void setValue( BaseValue currentValue ) {
        _existsValue = true;
        _value = currentValue;
    };

    friend class GenericMaterialProperty;
};

/** @typedef Definition d'une propriete materiau de type double */
typedef ElementaryMaterialProperty< ASTERDOUBLE > ElementaryMaterialPropertyReal;
/** @typedef Definition d'une propriete materiau de type double */
typedef ElementaryMaterialProperty< ASTERCOMPLEX > ElementaryMaterialPropertyComplex;
/** @typedef Definition d'une propriete materiau de type string */
typedef ElementaryMaterialProperty< std::string > ElementaryMaterialPropertyString;
/** @typedef Definition d'une propriete materiau de type Function */
typedef ElementaryMaterialProperty< FunctionPtr > ElementaryMaterialPropertyFunction;
/** @typedef Definition d'une propriete materiau de type Table */
typedef ElementaryMaterialProperty< TablePtr > ElementaryMaterialPropertyTable;
/** @typedef Definition d'une propriete materiau de type Function2D */
typedef ElementaryMaterialProperty< Function2DPtr > ElementaryMaterialPropertyFunction2D;
/** @typedef Definition d'une propriete materiau de type Formula */
typedef ElementaryMaterialProperty< FormulaPtr > ElementaryMaterialPropertyFormula;
/** @typedef Definition d'une propriete materiau de type DataStructure */
typedef ElementaryMaterialProperty< GenericFunctionPtr >
    ElementaryMaterialPropertyDataStructure;
/** @typedef Definition d'une propriete materiau de type vector double */
typedef ElementaryMaterialProperty< VectorReal > ElementaryMaterialPropertyVectorReal;
/** @typedef Definition d'une propriete materiau de type vector Function */
typedef ElementaryMaterialProperty< std::vector< FunctionPtr > >
    ElementaryMaterialPropertyVectorFunction;
/** @typedef Definition d'une propriete materiau de type Convertible string double */
typedef ElementaryMaterialProperty< StringToValueReal > ElementaryMaterialPropertyConvertible;

/**
 * @class GenericMaterialProperty
 * @brief Cette classe permet de definir un ensemble de type elementaire de propriete materielle
 * @author Nicolas Sellenet
 */
class GenericMaterialProperty {
  protected:
    /** @typedef std::map d'une chaine et d'un ElementaryMaterialPropertyReal */
    typedef std::map< std::string, ElementaryMaterialPropertyReal > mapStrEMPD;
    /** @typedef Iterateur sur mapStrEMPD */
    typedef mapStrEMPD::iterator mapStrEMPDIterator;
    typedef mapStrEMPD::const_iterator mapStrEMPDConstIterator;
    /** @typedef Valeur contenue dans un mapStrEMPD */
    typedef mapStrEMPD::value_type mapStrEMPDValue;

    /** @typedef std::map d'une chaine et d'un ElementaryMaterialPropertyComplex */
    typedef std::map< std::string, ElementaryMaterialPropertyComplex > mapStrEMPC;
    /** @typedef Iterateur sur mapStrEMPC */
    typedef mapStrEMPC::iterator mapStrEMPCIterator;
    typedef mapStrEMPC::const_iterator mapStrEMPCConstIterator;
    /** @typedef Valeur contenue dans un mapStrEMPC */
    typedef mapStrEMPC::value_type mapStrEMPCValue;

    /** @typedef std::map d'une chaine et d'un ElementaryMaterialPropertyString */
    typedef std::map< std::string, ElementaryMaterialPropertyString > mapStrEMPS;
    /** @typedef Iterateur sur mapStrEMPS */
    typedef mapStrEMPS::iterator mapStrEMPSIterator;
    typedef mapStrEMPS::const_iterator mapStrEMPSConstIterator;
    /** @typedef Valeur contenue dans un mapStrEMPS */
    typedef mapStrEMPS::value_type mapStrEMPSValue;

    /** @typedef std::map d'une chaine et d'un ElementaryMaterialPropertyDataStructure */
    typedef std::map< std::string, ElementaryMaterialPropertyDataStructure > mapStrEMPF;
    /** @typedef Iterateur sur mapStrEMPF */
    typedef mapStrEMPF::iterator mapStrEMPFIterator;
    typedef mapStrEMPF::const_iterator mapStrEMPFConstIterator;
    /** @typedef Valeur contenue dans un mapStrEMPF */
    typedef mapStrEMPF::value_type mapStrEMPFValue;

    /** @typedef std::map d'une chaine et d'un ElementaryMaterialPropertyTable */
    typedef std::map< std::string, ElementaryMaterialPropertyTable > mapStrEMPT;
    /** @typedef Iterateur sur mapStrEMPT */
    typedef mapStrEMPT::iterator mapStrEMPTIterator;
    typedef mapStrEMPT::const_iterator mapStrEMPTConstIterator;
    /** @typedef Valeur contenue dans un mapStrEMPT */
    typedef mapStrEMPT::value_type mapStrEMPTValue;

    /** @typedef std::map d'une chaine et d'un ElementaryMaterialPropertyVectorReal */
    typedef std::map< std::string, ElementaryMaterialPropertyVectorReal > mapStrEMPVD;
    /** @typedef Iterateur sur mapStrEMPVD */
    typedef mapStrEMPVD::iterator mapStrEMPVDIterator;
    typedef mapStrEMPVD::const_iterator mapStrEMPVDConstIterator;
    /** @typedef Valeur contenue dans un mapStrEMPVD */
    typedef mapStrEMPVD::value_type mapStrEMPVDValue;

    /** @typedef std::map d'une chaine et d'un ElementaryMaterialPropertyVectorFunction */
    typedef std::map< std::string, ElementaryMaterialPropertyVectorFunction > mapStrEMPVF;
    /** @typedef Iterateur sur mapStrEMPVF */
    typedef mapStrEMPVF::iterator mapStrEMPVFIterator;
    typedef mapStrEMPVF::const_iterator mapStrEMPVFConstIterator;
    /** @typedef Valeur contenue dans un mapStrEMPVF */
    typedef mapStrEMPVF::value_type mapStrEMPVFValue;

    /** @typedef std::map d'une chaine et d'un ElementaryMaterialPropertyConvertible */
    typedef std::map< std::string, ElementaryMaterialPropertyConvertible > mapStrEMPCSD;
    /** @typedef Iterateur sur mapStrEMPCSD */
    typedef mapStrEMPCSD::iterator mapStrEMPCSDIterator;
    typedef mapStrEMPCSD::const_iterator mapStrEMPCSDConstIterator;
    /** @typedef Valeur contenue dans un mapStrEMPCSD */
    typedef mapStrEMPCSD::value_type mapStrEMPCSDValue;

    /** @typedef std::list< std::string > */
    typedef std::list< std::string > ListString;
    typedef ListString::iterator ListStringIter;
    typedef ListString::const_iterator ListStringConstIter;

    /** @brief Chaine correspondant au nom Aster du MaterialProperty */
    // ex : ELAS ou ELASFo
    std::string _asterName;
    std::string _asterNewName;
    /** @brief Map contenant les noms des proprietes double ainsi que les
               MaterialProperty correspondant */
    mapStrEMPD _mapOfRealMaterialProperties;
    /** @brief Map contenant les noms des proprietes complex ainsi que les
               MaterialProperty correspondant */
    mapStrEMPC _mapOfComplexMaterialProperties;
    /** @brief Map contenant les noms des proprietes chaine ainsi que les
               MaterialProperty correspondant */
    mapStrEMPS _mapOfStringMaterialProperties;
    /** @brief Map contenant les noms des proprietes function ainsi que les
               MaterialProperty correspondant */
    mapStrEMPF _mapOfFunctionMaterialProperties;
    /** @brief Map contenant les noms des proprietes table ainsi que les
               MaterialProperty correspondant */
    mapStrEMPT _mapOfTableMaterialProperties;
    /** @brief Map contenant les noms des proprietes vector double ainsi que les
               MaterialProperty correspondant */
    mapStrEMPVD _mapOfVectorRealMaterialProperties;
    /** @brief Map contenant les noms des proprietes vector Function ainsi que les
               MaterialProperty correspondant */
    mapStrEMPVF _mapOfVectorFunctionMaterialProperties;
    /** @brief Map contenant les noms des proprietes double ainsi que les
               MaterialProperty correspondant */
    mapStrEMPCSD _mapOfConvertibleMaterialProperties;
    /** @brief Liste contenant les infos du .ORDR */
    VectorString _vectOrdr;
    /** @brief Vector of ordered keywords */
    VectorString _vectKW;
    /** @brief Has a traction function to build */
    bool _hasTractionFunction;


  public:
    /**
     * @brief Constructeur
     */
    GenericMaterialProperty() : _asterNewName( "" ), _hasTractionFunction(false) {};

    /**
     * @brief Constructeur
     */
    GenericMaterialProperty( const std::string asterName,
                                      const std::string asterNewName = "" )
        : _asterName( asterName ), _asterNewName( asterNewName ),
          _hasTractionFunction(false){};

    /**
     * @brief Recuperation du nom Aster du GenericMaterialProperty
     *        ex : 'ELAS', 'ELASFo', ...
     * @return Chaine contenant le nom Aster
     */
    const std::string getAsterName() const { return _asterName; };

    const std::string getAsterNewName() const { return _asterNewName; };

    /**
     * @brief Get number of properties containing a list of doubles
     */
    int getNumberOfListOfPropertiesReal() const {
        return _mapOfVectorRealMaterialProperties.size();
    };

    /**
     * @brief Get number of properties containing a list of functions
     */
    int getNumberOfListOfPropertiesFunction() const {
        return _mapOfVectorFunctionMaterialProperties.size();
    };

    ASTERDOUBLE getValueReal( std::string nameOfProperty )
    {
        auto curIter = _mapOfRealMaterialProperties.find( nameOfProperty );
        if ( curIter != _mapOfRealMaterialProperties.end() )
            return ( *curIter ).second.getValue();
        throw std::runtime_error( nameOfProperty + " is not a double value" );
    };

    ASTERCOMPLEX getValueComplex( std::string nameOfProperty )
    {
        auto curIter = _mapOfComplexMaterialProperties.find( nameOfProperty );
        if ( curIter != _mapOfComplexMaterialProperties.end() )
            return ( *curIter ).second.getValue();
        throw std::runtime_error( nameOfProperty + " is not a complex value" );
    };

    std::string getValueString( std::string nameOfProperty )
    {
        auto curIter = _mapOfStringMaterialProperties.find( nameOfProperty );
        if ( curIter != _mapOfStringMaterialProperties.end() )
            return ( *curIter ).second.getValue();
        throw std::runtime_error( nameOfProperty + " is not a string value" );
    };

    GenericFunctionPtr getValueGenericFunction( std::string nameOfProperty )
    {
        auto curIter = _mapOfFunctionMaterialProperties.find( nameOfProperty );
        if ( curIter != _mapOfFunctionMaterialProperties.end() )
            return ( *curIter ).second.getValue();
        throw std::runtime_error( nameOfProperty + " is not a function value" );
    };

    TablePtr getValueTable( std::string nameOfProperty )
    {
        auto curIter = _mapOfTableMaterialProperties.find( nameOfProperty );
        if ( curIter != _mapOfTableMaterialProperties.end() )
            return ( *curIter ).second.getValue();
        throw std::runtime_error( nameOfProperty + " is not a table value" );
    };

    bool hasValueReal( std::string nameOfProperty )
    {
        auto curIter = _mapOfRealMaterialProperties.find( nameOfProperty );
        if ( curIter == _mapOfRealMaterialProperties.end() )
            return false;
        return true;
    };

    bool hasValueComplex( std::string nameOfProperty )
    {
        auto curIter = _mapOfComplexMaterialProperties.find( nameOfProperty );
        if ( curIter == _mapOfComplexMaterialProperties.end() )
            return false;
        return true;
    };

    bool hasValueString( std::string nameOfProperty )
    {
        auto curIter = _mapOfStringMaterialProperties.find( nameOfProperty );
        if ( curIter == _mapOfStringMaterialProperties.end() )
            return false;
        return true;
    };

    bool hasValueGenericFunction( std::string nameOfProperty )
    {
        auto curIter = _mapOfFunctionMaterialProperties.find( nameOfProperty );
        if ( curIter == _mapOfFunctionMaterialProperties.end() )
            return false;
        return true;
    };

    bool hasValueTable( std::string nameOfProperty )
    {
        auto curIter = _mapOfTableMaterialProperties.find( nameOfProperty );
        if ( curIter == _mapOfTableMaterialProperties.end() )
            return false;
        return true;
    };

    /**
     * @brief Get number of properties with a value
     */
    int getNumberOfPropertiesWithValue() const;

    /**
     * @brief Fonction servant a fixer un parametre materiau au GenericMaterialProperty
     * @param nameOfProperty Nom de la propriete
     * @param value Real correspondant a la valeur donnee par l'utilisateur
     * @return Booleen valant true si la tache s'est bien deroulee
     */
    bool setValue( std::string nameOfProperty, ASTERDOUBLE value );

    /**
     * @brief Fonction servant a fixer un parametre materiau au GenericMaterialProperty
     * @param nameOfProperty Nom de la propriete
     * @param value Real correspondant a la valeur donnee par l'utilisateur
     * @return Booleen valant true si la tache s'est bien deroulee
     */
    bool setValue( std::string nameOfProperty, ASTERCOMPLEX value );

    /**
     * @brief Fonction servant a fixer un parametre materiau au GenericMaterialProperty
     * @param nameOfProperty Nom de la propriete
     * @param value string correspondant a la valeur donnee par l'utilisateur
     * @return Booleen valant true si la tache s'est bien deroulee
     */
    bool setValue( std::string nameOfProperty, std::string value );

    /**
     * @brief Fonction servant a fixer un parametre materiau au GenericMaterialProperty
     * @param nameOfProperty Nom de la propriete
     * @param value Function correspondant a la valeur donnee par l'utilisateur
     * @return Booleen valant true si la tache s'est bien deroulee
     */
    bool setValue( std::string nameOfProperty, FunctionPtr value );

    /**
     * @brief Fonction servant a fixer un parametre materiau au GenericMaterialProperty
     * @param nameOfProperty Nom de la propriete
     * @param value Table correspondant a la valeur donnee par l'utilisateur
     * @return Booleen valant true si la tache s'est bien deroulee
     */
    bool setValue( std::string nameOfProperty, TablePtr value );

    /**
     * @brief Fonction servant a fixer un parametre materiau au GenericMaterialProperty
     * @param nameOfProperty Nom de la propriete
     * @param value Function2D correspondant a la valeur donnee par l'utilisateur
     * @return Booleen valant true si la tache s'est bien deroulee
     */
    bool setValue( std::string nameOfProperty, Function2DPtr value );

    /**
     * @brief Fonction servant a fixer un parametre materiau au GenericMaterialProperty
     * @param nameOfProperty Nom de la propriete
     * @param value Formula correspondant a la valeur donnee par l'utilisateur
     * @return Booleen valant true si la tache s'est bien deroulee
     */
    bool setValue( std::string nameOfProperty, FormulaPtr value );

    /**
     * @brief Fonction servant a fixer un parametre materiau au GenericMaterialProperty
     * @param nameOfProperty Nom de la propriete
     * @param value VectorReal correspondant a la valeur donnee par l'utilisateur
     * @return Booleen valant true si la tache s'est bien deroulee
     */
    bool setValue( std::string nameOfProperty, VectorReal value );

    /**
     * @brief Fonction servant a fixer un parametre materiau au GenericMaterialProperty
     * @param nameOfProperty Nom de la propriete
     * @param value VectorFunction correspondant a la valeur donnee par l'utilisateur
     * @return Booleen valant true si la tache s'est bien deroulee
     */
    bool setValue( std::string nameOfProperty, VectorFunction value );

    /**
     * @brief Fonction permettant de définir le .ORDR
     * @param values Vecteur correspondant au mot clé ORDRE_PARAM
     * @return Booleen valant true si la tache s'est bien deroulee
     */
    bool setSortedListParameters( VectorString values ) {
        _vectOrdr = values;
        return true;
    };

    /**
     * @brief Construction du GenericMaterialProperty
     * @return Booleen valant true si la tache s'est bien deroulee
     * @todo vérifier les valeurs réelles par défaut du .VALR
     */
    virtual bool buildJeveuxVectors( JeveuxVectorComplex &complexValues,
                                     JeveuxVectorReal &doubleValues,
                                     JeveuxVectorChar16 &char16Values, JeveuxVectorChar16 &ordr,
                                     JeveuxVectorLong &kOrdr, std::vector< JeveuxVectorReal > &,
                                     std::vector< JeveuxVectorChar8 > & );

    /**
     * @brief Build ".RDEP" if necessary
     * @return true
     */
    virtual bool computeTractionFunction( FunctionPtr &doubleValues ) const;

    /**
     * @brief Function to know if ".RDEP" is necessary
     * @return true if ".RDEP" is necessary
     */
    bool hasTractionFunction() const { return _hasTractionFunction; };

    /**
     * @brief Function to know if ".RDEP" is necessary
     * @return true if ".RDEP" is necessary
     */
    void hasTractionFunction(const bool has_traction)
    { _hasTractionFunction = has_traction; };

    /**
     * @brief Function to know if behaviour own a list of double parameter
     */
    bool hasVectorOfRealParameters() const {
        if ( _mapOfVectorRealMaterialProperties.size() != 0 )
            return true;
        return false;
    };

    /**
     * @brief Function to know if behaviour own a list of double parameter
     */
    bool hasVectorOfFunctionParameters() const {
        if ( _mapOfVectorFunctionMaterialProperties.size() != 0 )
            return true;
        return false;
    };

  protected:
    bool addPropertyReal( std::string key, ElementaryMaterialPropertyReal value ) {
        _mapOfRealMaterialProperties[key] = value;
        _vectKW.push_back( key );
        _vectKW.push_back( value.getName() );
        return true;
    };

    bool addPropertyComplex( std::string key, ElementaryMaterialPropertyComplex value ) {
        _mapOfComplexMaterialProperties[key] = value;
        _vectKW.push_back( key );
        _vectKW.push_back( value.getName() );
        return true;
    };

    bool addPropertyString( std::string key, ElementaryMaterialPropertyString value ) {
        _mapOfStringMaterialProperties[key] = value;
        _vectKW.push_back( key );
        _vectKW.push_back( value.getName() );
        return true;
    };

    bool addPropertyFunction( std::string key, ElementaryMaterialPropertyDataStructure value ) {
        _mapOfFunctionMaterialProperties[key] = value;
        _vectKW.push_back( key );
        _vectKW.push_back( value.getName() );
        return true;
    };

    bool addPropertyTable( std::string key, ElementaryMaterialPropertyTable value ) {
        _mapOfTableMaterialProperties[key] = value;
        _vectKW.push_back( key );
        _vectKW.push_back( value.getName() );
        return true;
    };

    bool addPropertyVectorOfReal( std::string key,
                                    ElementaryMaterialPropertyVectorReal value ) {
        _mapOfVectorRealMaterialProperties[key] = value;
        _vectKW.push_back( key );
        _vectKW.push_back( value.getName() );
        return true;
    };

    bool addPropertyVectorOfFunction( std::string key,
                                      ElementaryMaterialPropertyVectorFunction value ) {
        _mapOfVectorFunctionMaterialProperties[key] = value;
        _vectKW.push_back( key );
        _vectKW.push_back( value.getName() );
        return true;
    };

    bool addConvertibleProperty( std::string key, ElementaryMaterialPropertyConvertible value ) {
        _mapOfConvertibleMaterialProperties[key] = value;
        _vectKW.push_back( key );
        _vectKW.push_back( value.getName() );
        return true;
    };
};

/** @typedef Pointeur intellignet vers un comportement materiau quelconque */
typedef boost::shared_ptr< GenericMaterialProperty > GenericMaterialPropertyPtr;

#endif
