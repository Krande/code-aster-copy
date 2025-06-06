#ifndef FUNCTION2D_H_
#define FUNCTION2D_H_

/**
 * @file Function2D.h
 * @brief Fichier entete de la classe Function2D
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

#include "Functions/GenericFunction.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "Supervis/ResultNaming.h"

/**
 * @class Function2D
 * @brief Cette classe correspond à une nappe
 */
class Function2D : public GenericFunction {
  private:
    // Vecteur Jeveux '.PARA'
    JeveuxVectorReal _parameters;
    // Vecteur Jeveux '.VALE'
    JeveuxCollectionReal _value;
    /** @brief Python attributes for t_nappe */
    py::object _t_nappe;

  public:
    /**
     * @typedef Function2DPtr
     * @brief Pointeur intelligent vers un Function2D
     */
    typedef std::shared_ptr< Function2D > Function2DPtr;

    /**
     * @brief Constructeur
     */
    Function2D() : Function2D( ResultNaming::getNewResultName() ) {};

    /**
     * @brief Constructeur
     */
    Function2D( const std::string name )
        : GenericFunction( name, "NAPPE", "NAPPE" ),
          _parameters( JeveuxVectorReal( getName() + ".PARA" ) ),
          _value( JeveuxCollectionReal( getName() + ".VALE" ) ),
          _t_nappe( py::none() ) {};

    /**
     * @brief Return the parameters values of the function
     */
    const JeveuxVectorReal getParameters() const { return _parameters; }

    /**
     * @brief Return the values of the function
     */
    const JeveuxCollectionReal getValues() const { return _value; }

    /**
     * @brief Get the result name
     * @return  name of the result
     */
    std::string getResultName();

    /**
     * @brief Return the maximum number of points of the functions
     */
    ASTERINTEGER maximumSize() const;

    /**
     * @brief Return the total number of points of the functions
     */
    ASTERINTEGER size() const;

    /**
     * @brief Getter for nappe property
     */
    const py::object &getTNappe() const { return _t_nappe; }

    /**
     * @brief Setter for nappe property
     */
    void setTNappe( py::object &t_nappe ) { _t_nappe = t_nappe; }
};

/**
 * @typedef Function2DPtr
 * @brief Pointeur intelligent vers un Function2D
 */
using Function2DPtr = std::shared_ptr< Function2D >;

#endif /* FUNCTION2D_H_ */
