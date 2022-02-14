#ifndef ELEMENTARYVECTOR_H_
#define ELEMENTARYVECTOR_H_

/**
 * @file ElementaryVector.h
 * @brief Definition of elementary vectors
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

#include "astercxx.h"

#include "DataFields/ElementaryTerm.h"
#include "DataFields/FieldOnCells.h"
#include "DataFields/FieldOnNodes.h"
#include "DataStructures/DataStructure.h"
#include "Loads/ListOfLoads.h"
#include "Loads/PhysicalQuantity.h"
#include "MemoryManager/JeveuxVector.h"
#include "Numbering/DOFNumbering.h"
#include "Supervis/ResultNaming.h"

/**
 * @class BaseElementaryVector
 * @brief Generic class for sd_vect_elem
 */
class BaseElementaryVector : public DataStructure {
  protected:
    /** @brief Objects for computation of elementary terms */
    JeveuxVectorChar24 _rerr;
    JeveuxVectorChar24 _relr;

    /** @brief Flag for empty datastructure */
    bool _isEmpty;

    /** @brief Vectors of RESUELEM */
    std::vector< ElementaryTermRealPtr > _realVector;

    /** @brief Liste de charges */
    ListOfLoadsPtr _listOfLoads;

    /** @brief Link betwwen load and vector */
    NamesMapChar24 _corichRept;

  public:
    /** @brief Constructor with a name */
    BaseElementaryVector( const std::string name, const std::string type = "VECT_ELEM" )
        : DataStructure( name, 19, type ),
          _rerr( JeveuxVectorChar24( getName() + ".RERR" ) ),
          _relr( JeveuxVectorChar24( getName() + ".RELR" ) ),
          _isEmpty( true ),
          _listOfLoads( new ListOfLoads() ),
          _corichRept( NamesMapChar24( "&&CORICH.REPT" ) ){};

    /** @brief Constructor with automatic name */
    BaseElementaryVector() : BaseElementaryVector( ResultNaming::getNewResultName() ){};

  public:
    /** @brief Add a load to elementary vector */
    template < typename... Args >
    void addLoad( const Args &...a ) {
        _listOfLoads->addLoad( a... );
    };

    /** @brief Set type of vector */
    void setType( const std::string newType ) { DataStructure::setType( newType ); };

    /**
     * @brief Assembly with dofNume and time (for load)
     * @param dofNume object DOFNumbering
     */
    FieldOnNodesRealPtr assembleWithLoadFunctions( const BaseDOFNumberingPtr &dofNume,
                                                   const ASTERDOUBLE &time );

    /**
     * @brief Assembly with dofNume and apply functions from loads
     * @param dofNume object DOFNumbering
     */
    FieldOnNodesRealPtr assembleWithLoadFunctions( const BaseDOFNumberingPtr &dofNume ) {
        return assembleWithLoadFunctions( dofNume, 0. );
    };

    /**
     * @brief Assembly with dofNume
     * @param dofNume object DOFNumbering
     */
    FieldOnNodesRealPtr assemble( const BaseDOFNumberingPtr dofNume ) const;

    /**
     * @brief Assembly with dofNume and mask on cells
     * @param dofNume object DOFNumbering
     * @param maskCell FieldOnCells to apply mask on each cell
     * @param maskInve flag to inverse mask
     */
    FieldOnNodesRealPtr assembleWithMask( const BaseDOFNumberingPtr &dofNume,
                                          const FieldOnCellsLongPtr &maskCell,
                                          const int &maskInve );

    /**
     * @brief Detect state of datastructure
     * @return true if empty datastructure
     */
    bool isEmpty() { return _isEmpty; };

    /**
     * @brief Set state of datastructure
     * @param bEmpty flag for state of datastructure
     */
    void isEmpty( bool bEmpty ) { _isEmpty = bEmpty; };

    /**
     * @brief Set list of loads
     * @param currentList list of loads
     */
    void setListOfLoads( const ListOfLoadsPtr &currentList ) { _listOfLoads = currentList; };

    /** @brief Function to build ElementaryTerm */
    bool build();

    friend class DiscreteComputation;
};

/** @typedef BaseElementaryVectorPtr */
typedef boost::shared_ptr< BaseElementaryVector > BaseElementaryVectorPtr;

/**
 * @class ElementaryVector
 * @brief Class for sd_vect_elem template
 */
template < typename ValueType, PhysicalQuantityEnum PhysicalQuantity >
class ElementaryVector : public BaseElementaryVector {
  public:
    /** @typedef ElementaryVectorPtr */
    typedef boost::shared_ptr< ElementaryVector< ValueType, PhysicalQuantity > >
        ElementaryVectorPtr;

    /** @brief Constructor with a name */
    ElementaryVector( const std::string name )
        : BaseElementaryVector(
              name, "VECT_ELEM_" + std::string( PhysicalQuantityNames[PhysicalQuantity] ) +
                        ( typeid( ValueType ) == typeid( ASTERDOUBLE ) ? "_R" : "_C" ) ){};

    /** @brief Constructor with automatic name */
    ElementaryVector() : ElementaryVector( ResultNaming::getNewResultName() ){};
};

/** @typedef Elementary vector for displacement-double */
template class ElementaryVector< ASTERDOUBLE, Displacement >;
typedef ElementaryVector< ASTERDOUBLE, Displacement > ElementaryVectorDisplacementReal;
typedef boost::shared_ptr< ElementaryVectorDisplacementReal > ElementaryVectorDisplacementRealPtr;

/** @typedef Elementary vector for temperature-double */
template class ElementaryVector< ASTERDOUBLE, Temperature >;
typedef ElementaryVector< ASTERDOUBLE, Temperature > ElementaryVectorTemperatureReal;
typedef boost::shared_ptr< ElementaryVectorTemperatureReal > ElementaryVectorTemperatureRealPtr;

/** @typedef Elementary vector for pressure-complex */
template class ElementaryVector< ASTERCOMPLEX, Pressure >;
typedef ElementaryVector< ASTERCOMPLEX, Pressure > ElementaryVectorPressureComplex;
typedef boost::shared_ptr< ElementaryVectorPressureComplex > ElementaryVectorPressureComplexPtr;

#endif /* ELEMENTARYVECTOR_H_ */
