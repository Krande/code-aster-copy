#ifndef FUNCTION_H_
#define FUNCTION_H_

/**
 * @file Function.h
 * @brief Implementation of functions.
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

 * person_in_charge: mathieu.courtois@edf.fr
 */
#include "Functions/GenericFunction.h"
#include "MemoryManager/JeveuxVector.h"

/**
 * class BaseFunction
 *   Create a datastructure for a function with real values
 * @author Mathieu Courtois
 */
class BaseFunction : public GenericFunction {
  private:
  protected:
    // Vecteur Jeveux '.VALE'
    JeveuxVectorReal _value;

  public:
    /**
     * @typedef BaseFunctionPtr
     * @brief Pointeur intelligent vers un BaseFunction
     */
    typedef std::shared_ptr< BaseFunction > BaseFunctionPtr;

    /**
     * Constructeur
     */
    BaseFunction( const std::string type, const std::string type2 );

    BaseFunction( const std::string name, const std::string type, const std::string type2 );

    ~BaseFunction() {};

    /**
     * @brief Allocate function
     */
    virtual void allocate( ASTERINTEGER size );

    /**
     * @brief Deallocate function
     */
    void deallocate();

    /**
     * @brief Get the result name
     * @return  name of the result
     */
    std::string getResultName();

    /**
     * @brief Set the function type to CONSTANT
     */
    void setAsConstant();

    /**
     * @brief Definition of the name of the parameter (abscissa)
     * @param name name of the parameter
     * @type  name string
     */
    void setParameterName( const std::string name );

    /**
     * @brief Definition of the name of the result (ordinate)
     * @param name name of the result
     * @type  name string
     */
    void setResultName( const std::string name );

    /**
     * @brief Definition of the type of interpolation
     * @param interpolation type of interpolation
     * @type  interpolation string
     * @todo checking
     */
    void setInterpolation( const std::string type );

    /**
     * @brief Assign the values of the function
     * @param absc values of the abscissa
     * @type  absc vector of double
     * @param ord values of the ordinates
     * @type  ord vector of double
     */
    virtual void setValues( const VectorReal &absc, const VectorReal &ord );

    /**
     * @brief Return the values of the function
     */
    const JeveuxVectorReal getValues() const { return _value; }

    /**
     * @brief Return a pointer to the vector of data
     */
    const ASTERDOUBLE *getDataPtr() const {
        _value->updateValuePointer();
        return _value->getDataPtr();
    }

    /**
     * @brief Return the number of points of the function
     */
    virtual ASTERINTEGER maximumSize() const {
        if ( !_value.exists() )
            return 0;
        _value->updateValuePointer();
        return _value->size() / 2;
    }

    /**
     * @brief Return the number of points of the function
     */
    virtual ASTERINTEGER size() const { return maximumSize(); }

    /**
     * @brief Update the pointers to the Jeveux objects
     * @return Return true if ok
     */
    bool build() {
        _property->updateValuePointer();
        _value->updateValuePointer();
        return true;
    }
};

/**
 * class Function
 *   Create a datastructure for a function with real values
 * @author Mathieu Courtois
 */
class Function : public BaseFunction {

  public:
    /**
     * @typedef FunctionPtr
     * @brief Pointeur intelligent vers un Function
     */
    typedef std::shared_ptr< Function > FunctionPtr;

    /**
     * Constructeur
     */
    Function() : BaseFunction( "FONCTION", "FONCTION" ) {};

    Function( const std::string name ) : BaseFunction( name, "FONCTION", "FONCTION" ) {};
};

/**
 * class FunctionComplex
 *   Create a datastructure for a function with complex values
 * @author Mathieu Courtois
 */
class FunctionComplex : public BaseFunction {

  public:
    /**
     * @typedef FunctionPtr
     * @brief Pointeur intelligent vers un FunctionComplex
     */
    typedef std::shared_ptr< FunctionComplex > FunctionComplexPtr;

    /**
     * Constructeur
     */
    FunctionComplex( const std::string name ) : BaseFunction( name, "FONCTION_C", "FONCT_C" ) {};

    FunctionComplex() : BaseFunction( "FONCTION_C", "FONCT_C" ) {};

    /**
     * @brief Allocate function
     */
    void allocate( ASTERINTEGER size );

    /**
     * @brief Return the number of points of the function
     */
    virtual ASTERINTEGER maximumSize() const { return _value->size() / 3; }

    /**
     * @brief Return the number of points of the function
     */
    ASTERINTEGER size() const { return _value->size() / 3; }

    /**
     * @brief Assign the values of the function
     * @param absc values of the abscissa
     * @type  absc vector of double
     * @param ord values of the ordinates (real1, imag1, real2, imag2...)
     * @type  ord vector of double
     */
    void setValues( const VectorReal &absc, const VectorReal &ord );

    /**
     * @brief Assign the values of the function
     * @param absc values of the abscissa
     * @type  absc vector of double
     * @param ord values of the ordinates
     * @type  ord vector of complex
     */
    void setValues( const VectorReal &absc, const VectorComplex &ord );
};

/**
 * @typedef BaseFunctionPtr
 * @brief  Pointer to a BaseFunction
 */
typedef std::shared_ptr< BaseFunction > BaseFunctionPtr;

/**
 * @typedef FunctionPtr
 * @brief  Pointer to a Function
 */
typedef std::shared_ptr< Function > FunctionPtr;

/**
 * @typedef FunctionComplexPtr
 * @brief  Pointer to a FunctionComplex
 */
typedef std::shared_ptr< FunctionComplex > FunctionComplexPtr;

/**
 * @name emptyRealFunction
 * @brief  Empty function
 */
extern FunctionPtr emptyRealFunction;

/**
 * @typedef ListOfFunctions
 * @brief List of double functions
 */
typedef std::list< FunctionPtr > ListOfFunctions;

#endif /* FUNCTION_H_ */
