#ifndef ELEMENTARYCOMPUTE_H_
#define ELEMENTARYCOMPUTE_H_

/**
 * @file ElementaryCompute.h
 * @brief Header of class ElementaryCompute
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

#include "astercxx.h"

#include "Discretization/ElementaryCharacteristics.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"

/**
 * @class ElementaryCompute
 * @brief Defining objects to compute vect_elem and matr_elem
 */
class ElementaryCompute {
  private:
    /** @brief Objet Jeveux '.RERR' : descriptor*/
    JeveuxVectorChar24 _rerr;

    /** @brief Objet Jeveux '.RELR' : list of name of elementary results */
    JeveuxVectorChar24 _relr;

    /** @brief Option to compute */
    std::string _option;

  public:
    /** @brief Constructor by predefined name */
    ElementaryCompute( const std::string elemName, const std::string option ) : _option( option ) {
        std ::string baseName( elemName, 0, 19 );
        _rerr = JeveuxVectorChar24( baseName + ".RERR" );
        _relr = JeveuxVectorChar24( baseName + ".RELR" );
    };

    /** @brief Constructor by predefined name with calls to Fortran */
    ElementaryCompute( const std::string baseName )
        : ElementaryCompute( baseName, "WRAP_FORTRAN" ){};

    /** @brief Create descriptor */
    void createDescriptor( const ModelPtr &currentModel, const MaterialFieldPtr &currMaterialField,
                           const ElementaryCharacteristicsPtr &currElemChara );

    /** @brief Get option */
    std::string getOption() const { return _option; }

    void setOption( const std::string &option ) { _option = option; };

    /** @brief Get list of elementary terms */
    std::vector< JeveuxChar24 > getNameOfElementaryTerms() { return _relr->toVector(); };

    /** @brief Get descriptor */
    std::vector< JeveuxChar24 > getDescriptor() { return _rerr->toVector(); };

    /** @brief Add elementary term */
    void addElementaryTerm( const std::string elemTermName ) {
        if ( !hasElementaryTerm() ) {
            _relr->reserve( 10 );
        }

        _relr->updateValuePointer();
        _relr->push_back( elemTermName );
    };

    /** @brief Has elementary term ? */
    bool hasElementaryTerm() { return _relr->exists(); };

    /** @brief Get number of elementary term  */
    ASTERINTEGER getNumberOfElementaryTerms() const { return _relr->size(); };
};

using ElementaryComputePtr = std::shared_ptr< ElementaryCompute >;

#endif /* ELEMENTARYCOMPUTE_H_ */
