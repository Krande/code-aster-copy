#ifndef ELEMENTARYMATRIX_H_
#define ELEMENTARYMATRIX_H_

/**
 * @file ElementaryMatrix.h
 * @brief Definition of elementary matrices
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

#include "aster_fort_calcul.h"
#include "astercxx.h"

#include "DataFields/ElementaryTerm.h"
#include "LinearAlgebra/BaseElementaryMatrix.h"

/**
 * @class ElementaryMatrix
 * @brief Class for sd_matr_elem template
 */
template < typename ValueType, PhysicalQuantityEnum PhysicalQuantity >
class ElementaryMatrix : public BaseElementaryMatrix {
  private:
    /** @brief Vectors of RESUELEM */
    std::vector< boost::shared_ptr< ElementaryTerm< ValueType > > > _elemTerm;

  public:
    /** @typedef ElementaryMatrixPtr */
    typedef boost::shared_ptr< ElementaryMatrix< ValueType, PhysicalQuantity > >
        ElementaryMatrixPtr;

    /** @brief Constructor with a name */
    ElementaryMatrix( const std::string name )
        : BaseElementaryMatrix(
              name, "MATR_ELEM_" + std::string( PhysicalQuantityNames[PhysicalQuantity] ) +
                        ( typeid( ValueType ) == typeid( ASTERDOUBLE ) ? "_R" : "_C" ) ) {
        _elemTerm.clear();
    };

    /** @brief Constructor with automatic name */
    ElementaryMatrix() : ElementaryMatrix( ResultNaming::getNewResultName() ){};

    /**
     * @brief Function to update ElementaryTerm
     */
    bool build() {
        if ( _elemComp->hasElementaryTerm() ) {
            CALL_REDETR( getName() );

            SetString elemSave;
            for ( auto &elemTerm : _elemTerm ) {
                elemSave.insert( trim( elemTerm->getName() ) );
            }

            auto elemTermNames = _elemComp->getNameOfElementaryTerms();

            SetString elemKeep;
            for ( auto &elemTerm : elemTermNames ) {
                const std::string name = trim( elemTerm.toString() );
                elemKeep.insert( name );
                if ( name != " " && elemSave.count( name ) == 0 ) {
                    _elemTerm.push_back(
                        boost::make_shared< ElementaryTerm< ValueType > >( name ) );
                }
            }

            // clean ElementaryTerm
            std::vector< boost::shared_ptr< ElementaryTerm< ValueType > > > elemTermNew;
            elemTermNew.reserve( _elemTerm.size() );
            for ( auto &elemTerm : _elemTerm ) {
                auto name = trim( elemTerm->getName() );
                if ( elemKeep.count( name ) > 0 ) {
                    elemTermNew.push_back( elemTerm );
                }
            }
            _elemTerm = std::move( elemTermNew );

            if ( _elemTerm.size() > 0 && !isMPIFull() && !getMesh()->isParallel() ) {
                std::string type = "MATR_ELEM";
                CALLO_SDMPIC( type, getName() );
            }
        }
        _isEmpty = false;
        return true;
    };

    /**
     * @brief is MPI_COMPLET ?
     */
    bool isMPIFull() {
        for ( auto &elemTerm : _elemTerm ) {
            if ( !elemTerm->isMPIFull() ) {
                return false;
            }
        }

        return true;
    };

    /**
     * @brief Add elementary term
     */
    void addElementaryTerm( const boost::shared_ptr< ElementaryTerm< ValueType > > &elemTerm ) {
        _elemComp->addElementaryTerm( elemTerm->getName() );
        _elemTerm.push_back( elemTerm );
    };

    std::vector< boost::shared_ptr< ElementaryTerm< ValueType > > > getElementaryTerms() const {
        return _elemTerm;
    }

    /**
     * @brief Get all FiniteElementDescriptors
     * @return vector of all FiniteElementDescriptors
     */
    std::vector< FiniteElementDescriptorPtr > getFiniteElementDescriptors() {
        std::vector< FiniteElementDescriptorPtr > FEDs;
        FEDs.reserve( _elemTerm.size() );

        for ( auto &elemTerm : _elemTerm ) {
            auto FED = elemTerm->getFiniteElementDescriptor();
            if ( FED ) {
                FEDs.push_back( FED );
            }
        }

        return FEDs;
    };
};

/** @typedef Elementary matrix for displacement-double */
template class ElementaryMatrix< ASTERDOUBLE, Displacement >;
typedef ElementaryMatrix< ASTERDOUBLE, Displacement > ElementaryMatrixDisplacementReal;
typedef boost::shared_ptr< ElementaryMatrixDisplacementReal > ElementaryMatrixDisplacementRealPtr;

/** @typedef Elementary matrix for displacement-complex */
template class ElementaryMatrix< ASTERCOMPLEX, Displacement >;
typedef ElementaryMatrix< ASTERCOMPLEX, Displacement > ElementaryMatrixDisplacementComplex;
typedef boost::shared_ptr< ElementaryMatrixDisplacementComplex >
    ElementaryMatrixDisplacementComplexPtr;

/** @typedef Elementary matrix for temperature-double */
template class ElementaryMatrix< ASTERDOUBLE, Temperature >;
typedef ElementaryMatrix< ASTERDOUBLE, Temperature > ElementaryMatrixTemperatureReal;
typedef boost::shared_ptr< ElementaryMatrixTemperatureReal > ElementaryMatrixTemperatureRealPtr;

/** @typedef Elementary matrix for pressure-complex */
template class ElementaryMatrix< ASTERCOMPLEX, Pressure >;
typedef ElementaryMatrix< ASTERCOMPLEX, Pressure > ElementaryMatrixPressureComplex;
typedef boost::shared_ptr< ElementaryMatrixPressureComplex > ElementaryMatrixPressureComplexPtr;

#endif /* ELEMENTARYMATRIX_H_ */
