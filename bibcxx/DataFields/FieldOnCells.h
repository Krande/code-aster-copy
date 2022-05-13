#ifndef FIELDONCELLS_H_
#define FIELDONCELLS_H_

/**
 * @file FieldOnCells.h
 * @brief Header of class for FieldOnCells
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

#include "aster_fort_ds.h"
#include "aster_fort_superv.h"
#include "astercxx.h"

#include "Behaviours/BehaviourProperty.h"
#include "DataFields/DataField.h"
#include "DataFields/SimpleFieldOnCells.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "PythonBindings/LogicalUnitManager.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/Exceptions.h"

/** @brief Forward declaration of ElementaryCharacteristics */
class ElementaryCharacteristics;
using ElementaryCharacteristicsPtr = std::shared_ptr< ElementaryCharacteristics >;

/**
 * @class FieldOnCells
 * @brief Template class for FieldOnCells
 */
template < class ValueType >
class FieldOnCells : public DataField {
  private:
    using SimpleFieldOnCellsValueType = SimpleFieldOnCells< ValueType >;
    using SimpleFieldOnCellsValueTypePtr = std::shared_ptr< SimpleFieldOnCellsValueType >;

    /** @brief Vecteur Jeveux '.CELD' */
    JeveuxVectorLong _descriptor;
    /** @brief Vecteur Jeveux '.CELK' */
    JeveuxVectorChar24 _reference;
    /** @brief Vecteur Jeveux '.CELV' */
    JeveuxVector< ValueType > _valuesList;
    /** @brief Finite element description */
    FiniteElementDescriptorPtr _dofDescription;
    /** @brief Object for dynamic fields  (as VARI_ELGA) */
    SimpleFieldOnCellsLongPtr _DCEL;

  public:
    using FieldOnCellsPtr = std::shared_ptr< FieldOnCells >;

    /**
     * @brief Constructor
     * @param name Jeveux name of the field
     */
    FieldOnCells( const std::string name )
        : DataField( name, "CHAM_ELEM" ),
          _descriptor( JeveuxVectorLong( getName() + ".CELD" ) ),
          _reference( JeveuxVectorChar24( getName() + ".CELK" ) ),
          _valuesList( JeveuxVector< ValueType >( getName() + ".CELV" ) ),
          _dofDescription( nullptr ){};

    /** @brief Constructor with automatic name */
    FieldOnCells() : FieldOnCells( ResultNaming::getNewResultName() ){};

    /** @brief Constructor with automatic name and FE Descriptor*/
    FieldOnCells( const FiniteElementDescriptorPtr FEDesc )
        : FieldOnCells( ResultNaming::getNewResultName() ) {
        _dofDescription = FEDesc;
    };

    /**
     * @brief Constructor for empty FieldOnCells with dynamic components
     * @param model model
     * @param behaviour Description of behaviour (for size of dynamic components as VARI_ELGA)
     * @param carael Description of elementary characteristics (for size of dynamic components as
     * VARI_ELGA)
     * @param typcham Type de champ à calculer
     */
    FieldOnCells( const ModelPtr &model, const BehaviourPropertyPtr behaviour,
                  const std::string &typcham, const ElementaryCharacteristicsPtr carael = nullptr,
                  const FiniteElementDescriptorPtr FEDesc = nullptr );

    /**
     * @brief Constructor for empty FieldOnCells based on specific physical quantity
     */
    FieldOnCells( const ModelPtr &model, const std::string option, const std::string paraName,
                  const FiniteElementDescriptorPtr FEDesc = nullptr );

    /** @brief Copy constructor */
    FieldOnCells( const std::string &name, const FieldOnCells &toCopy ) : FieldOnCells( name ) {
        // JeveuxVector to be duplicated
        *( _descriptor ) = *( toCopy._descriptor );
        *( _reference ) = *( toCopy._reference );
        *( _valuesList ) = *( toCopy._valuesList );
        *( _title ) = *( toCopy._title );
        // Pointers to be copied
        setDescription( toCopy._dofDescription );
        updateValuePointers();
    }

    /** @brief Move constructor */
    FieldOnCells( FieldOnCells &&other ) : DataField{ std::move( other ) } {
        _descriptor = other._descriptor;
        _reference = other._reference;
        _valuesList = other._valuesList;
        _title = other._title;
        setDescription( other._dofDescription );
        updateValuePointers();
    }

    /**
     * @brief Copy constructor
     */
    FieldOnCells( const FieldOnCells &toCopy )
        : FieldOnCells( DataStructureNaming::getNewName(), toCopy ){};

    /**
     * @brief Wrap of copy constructor
     */
    FieldOnCells duplicate() { return *this; }

    /**
     * @brief
     * @return
     */
    void deallocate() {
        _descriptor->deallocate();
        _reference->deallocate();
        _valuesList->deallocate();
        _dofDescription = nullptr;
    };

    /**
     * @brief
     * @return
     */
    SimpleFieldOnCellsValueTypePtr exportToSimpleFieldOnCells() {
        SimpleFieldOnCellsValueTypePtr toReturn( new SimpleFieldOnCellsValueType() );
        const std::string resultName = toReturn->getName();
        const std::string inName = getName();
        const std::string copyNan( "OUI" );
        CALLO_CELCES_WRAP( inName, JeveuxMemoryTypesNames[Permanent], resultName );
        toReturn->updateValuePointers();
        return toReturn;
    };

    /**
     * @brief Check if fields are OK for +, +=, ...
     * @return true if compatible
     */
    ASTERBOOL isSimilarTo( const FieldOnCells< ValueType > &tmp2 ) const {
        bool similar = ( _descriptor->size() == tmp2._descriptor->size() );
        similar = ( similar && ( _reference->size() == tmp2._reference->size() ) );
        similar = ( similar && ( _valuesList->size() == tmp2._valuesList->size() ) );
        return similar;
    }

    /**
     * @brief Get the mesh
     */
    BaseMeshPtr getMesh() const {
        if ( _dofDescription ) {
            return _dofDescription->getMesh();
        }

        return nullptr;
    };

    ModelPtr getModel() const {
        if ( _dofDescription ) {
            return _dofDescription->getModel();
        }

        return nullptr;
    }

    /**
     * @brief Set the description of finite elements
     * @param curDesc object FiniteElementDescriptorPtr
     */
    void setDescription( const FiniteElementDescriptorPtr &curDesc ) {
        if ( !curDesc && _dofDescription ) {
            AS_ABORT( "FiniteElementDescriptor is empty" );
        }

        if ( _dofDescription && curDesc && _dofDescription != curDesc ) {
            std::string msg =
                "FiniteElementDescriptor inconsistents: " + _dofDescription->getName() + " vs " +
                curDesc->getName();
            AS_ABORT( msg );
        }

        _dofDescription = curDesc;
    };

    /**
     * @brief Get the description of finite elements
     */
    FiniteElementDescriptorPtr getDescription() const { return _dofDescription; };

    /**
     * @brief Update field and build FiniteElementDescriptor if necessary
     */
    ASTERBOOL build() {
        if ( !_dofDescription ) {
            CALL_JEMARQ();
            _reference->updateValuePointer();
            const std::string ligrel = trim( ( *_reference )[0].toString() );
            CALL_JEDEMA();

            if ( ligrel.substr( 0, 8 ) == getName().substr( 0, 8 ) ) {
                setDescription( std::make_shared< FiniteElementDescriptor >( ligrel, getMesh() ) );
            }
        }
        return true;
    };

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
     */
    void updateValuePointers() const {
        _descriptor->updateValuePointer();
        _reference->updateValuePointer();
        _valuesList->updateValuePointer();
    };

    /**
     * @brief Transormer les valeurs de _valuesList en appliquant
     *         la fonction "func" à chaque valeur
     * @return renvoie un nouveau objet de FieldOnCells
     *         avec les valeurs transformées
     */
    template < class type = ValueType >
    typename std::enable_if< std::is_same< type, ASTERDOUBLE >::value,
                             FieldOnCells< ValueType > >::type
    transform( py::object &func ) {
        if ( !PyCallable_Check( func.ptr() ) )
            raiseAsterError( "Input parameter to the transform \
        method should be a callable Python object" );

        FieldOnCells< ValueType > tmp( *this );
        updateValuePointers();

        ASTERINTEGER size = _valuesList->size();
        for ( auto i = 0; i < size; i++ ) {
            PyObject *res = PyObject_CallFunction( func.ptr(), "d", ( *_valuesList )[i] );
            if ( PyFloat_Check( res ) ) {
                tmp[i] = (ASTERDOUBLE)PyFloat_AsDouble( res );
            } else {
                PyErr_Format( PyExc_ValueError, "Returned value of \
                    type different from ASTERDOUBLE" );
                PyErr_Print();
            }
            Py_XDECREF( res );
        }

        return tmp;
    };

    template < class type = ValueType >
    typename std::enable_if< std::is_same< type, ASTERCOMPLEX >::value,
                             FieldOnCells< ValueType > >::type
    transform( py::object &func ) {
        if ( !PyCallable_Check( func.ptr() ) )
            raiseAsterError( "Input parameter to the transform \
        method should be a callable Python object" );

        FieldOnCells< ValueType > tmp( *this );
        _valuesList->updateValuePointer();

        ASTERINTEGER size = _valuesList->size();

        Py_complex val;
        for ( auto i = 0; i < size; i++ ) {
            val.real = ( *_valuesList )[i].real();
            val.imag = ( *_valuesList )[i].imag();
            PyObject *res = PyObject_CallFunction( func.ptr(), "D", val );
            if ( PyComplex_Check( res ) ) {
                ASTERDOUBLE re = (ASTERDOUBLE)PyComplex_RealAsDouble( res );
                ASTERDOUBLE im = (ASTERDOUBLE)PyComplex_ImagAsDouble( res );
                tmp[i] = { re, im };
            } else {
                PyErr_Format( PyExc_ValueError, "Returned value of \
                    type different from ASTERCOMPLEX" );
                PyErr_Print();
            }
            // Py_DECREF(res);
            Py_XDECREF( res );
        }
        return tmp;
    };

    // OVERLOADING C++ OPERATORS

    /**
     * @brief Unary Minus overloading
     * @return Updated field
     */
    FieldOnCells< ValueType > operator-() const {
        FieldOnCells< ValueType > tmp( *this );
        ( *tmp._valuesList ) *= ValueType( -1 );
        return tmp;
    };

    /**
     * @brief Shorthand + operator assignement
     * @return Updated field
     */
    FieldOnCells< ValueType > &operator+=( const FieldOnCells< ValueType > &rhs ) {
        if ( !this->isSimilarTo( rhs ) )
            raiseAsterError( "Fields have incompatible shapes" );
        ( *_valuesList ) += ( *rhs._valuesList );
        return *this;
    };

    /**
     * @brief Shorthand - operator assignement
     * @return Updated field
     */
    FieldOnCells< ValueType > &operator-=( const FieldOnCells< ValueType > &rhs ) {
        if ( !this->isSimilarTo( rhs ) )
            raiseAsterError( "Fields have incompatible shapes" );
        ( *_valuesList ) -= ( *rhs._valuesList );
        return *this;
    };

    /**
     * @brief Subscripting overloading
     * @param i subscript
     * @return value at position i
     */
    ValueType &operator[]( int i ) { return _valuesList->operator[]( i ); };

    const ValueType &operator[]( int i ) const {
        return const_cast< ValueType & >(
            const_cast< FieldOnCells< ValueType > * >( this )->operator[]( i ) );
    };

    /**
     * @brief Plus overloading
     * @return New field
     */
    friend FieldOnCells< ValueType > operator+( FieldOnCells< ValueType > lhs,
                                                const FieldOnCells< ValueType > &rhs ) {
        if ( !lhs.isSimilarTo( rhs ) )
            raiseAsterError( "Fields have incompatible shapes" );
        lhs += rhs;
        return lhs;
    };

    /**
     * @brief Minus overloading
     * @return New field
     */
    friend FieldOnCells< ValueType > operator-( FieldOnCells< ValueType > lhs,
                                                const FieldOnCells< ValueType > &rhs ) {
        if ( !lhs.isSimilarTo( rhs ) )
            raiseAsterError( "Fields have incompatible shapes" );
        lhs -= rhs;
        return lhs;
    };

    /**
     * @brief Multiply by a scalar on right overloading
     * @return New field
     */

    friend FieldOnCells< ValueType > operator*( FieldOnCells< ValueType > lhs,
                                                const ASTERDOUBLE &scal ) {
        ( *lhs._valuesList ) *= scal;
        return lhs;
    };

    /**
     * @brief Multiply by a scalar on left overloading
     * @return New field
     */

    friend FieldOnCells< ValueType > operator*( const ASTERDOUBLE &scal,
                                                const FieldOnCells< ValueType > &rhs ) {
        return rhs * scal;
    };

    // some getters

    /**
     * @brief Get values of the field
     *
     */
    const JeveuxVector< ValueType > &getValues() const { return _valuesList; }

    /**
     * @brief Set the Values object
     *
     * @param value Value to affect
     */
    void setValues( const ValueType &value ) {
        _valuesList->updateValuePointer();
        _valuesList->assign( value );
    };

    /**
     * @brief Size of the FieldOnNodes
     */
    const ASTERINTEGER size() const { return _valuesList->size(); }

    // norm and dot methods

    /**
     * @brief Comput norm
     * @param normType Type of norm ("NORM_1","NORM_2","NORM_INFINITY")
     */
    template < class type = ValueType >
    typename std::enable_if< std::is_same< type, ASTERDOUBLE >::value, ASTERDOUBLE >::type
    norm( const std::string normType ) const {
        ASTERDOUBLE norme = 0.0;
        ASTERINTEGER beg = 0, end = 0, nbgrp = 0;

        const int rank = getMPIRank();

        _valuesList->updateValuePointer();
        _descriptor->updateValuePointer();

        JeveuxVectorLong CellsRank = getMesh()->getCellsRank();
        CellsRank->updateValuePointer();

        JeveuxCollectionLong collec = _dofDescription->getListOfGroupOfCells();
        JeveuxVectorLong descr = _descriptor;
        nbgrp = ( *descr )[1];

        for ( auto i = 0; i < nbgrp; i++ ) {

            ASTERINTEGER adress = ( *descr )[4 + i];
            if ( ( *descr )[adress + 2] == 0 )
                continue;

            ASTERINTEGER nel = ( *descr )[adress];
            auto liel = ( *collec )[i + 1];

            if ( normType == "NORM_1" ) {
                for ( auto p = 0; p < nel; p++ ) {

                    if ( ( *CellsRank )[liel[p] - 1] != rank )
                        continue;
                    beg = ( *descr )[adress + 3 + 4 * p + 4] - 1;
                    end = beg + ( *descr )[adress + 3 + 4 * p + 3];

                    for ( int pos = beg; pos < end; ++pos ) {
                        norme += std::abs( ( *this )[pos] );
                    }
                }
            } else if ( normType == "NORM_2" ) {
                for ( auto p = 0; p < nel; p++ ) {

                    if ( ( *CellsRank )[liel[p] - 1] != rank )
                        continue;
                    beg = ( *descr )[adress + 3 + 4 * p + 4] - 1;
                    end = beg + ( *descr )[adress + 3 + 4 * p + 3];

                    for ( int pos = beg; pos < end; ++pos ) {
                        norme += ( *this )[pos] * ( *this )[pos];
                    }
                }
            } else if ( normType == "NORM_INFINITY" ) {
                for ( auto p = 0; p < nel; p++ ) {

                    if ( ( *CellsRank )[liel[p] - 1] != rank )
                        continue;
                    beg = ( *descr )[adress + 3 + 4 * p + 4] - 1;
                    end = beg + ( *descr )[adress + 3 + 4 * p + 3];

                    for ( int pos = beg; pos < end; ++pos ) {
                        norme = std::max( norme, std::abs( ( *this )[pos] ) );
                    }
                }
            } else {
                AS_ASSERT( false );
            }
        }

#ifdef ASTER_HAVE_MPI
        if ( getMesh()->isParallel() ) {
            ASTERDOUBLE norm2 = norme;
            if ( normType == "NORM_1" || normType == "NORM_2" )
                AsterMPI::all_reduce( norm2, norme, MPI_SUM );
            else
                AsterMPI::all_reduce( norm2, norme, MPI_MAX );
        }
#endif

        if ( normType == "NORM_2" )
            norme = std::sqrt( norme );

        return norme;
    }

    /**
     * @brief Dot product
     * @param tmp object FieldOnCellsPtr
     */
    template < class type = ValueType >
    typename std::enable_if< std::is_same< type, ASTERDOUBLE >::value, ASTERDOUBLE >::type
    dot( const FieldOnCellsPtr &tmp ) const {
        tmp->updateValuePointers();
        _valuesList->updateValuePointer();
        ASTERINTEGER taille = _valuesList->size();

        if ( taille != tmp->size() )
            raiseAsterError( "Incompatible size" );

        ASTERDOUBLE ret = 0.0;
        for ( auto pos = 0; pos < taille; ++pos ) {
            ret += ( *this )[pos] * ( *tmp )[pos];
        }
        return ret;
    }

    bool printMedFile( const std::string fileName, bool local = true ) const;

    friend class FieldBuilder;
};

template < class ValueType >
bool FieldOnCells< ValueType >::printMedFile( const std::string fileName, bool local ) const {
    LogicalUnitFile a( fileName, Binary, New );
    int retour = a.getLogicalUnit();
    CommandSyntax cmdSt( "IMPR_RESU" );

    SyntaxMapContainer dict;
    dict.container["FORMAT"] = "MED";
    dict.container["UNITE"] = (ASTERINTEGER)retour;

    if ( getMesh()->isParallel() ) {
        dict.container["PROC0"] = "NON";
        if ( !local )
            dict.container["FICHIER_UNIQUE"] = "OUI";
    } else
        dict.container["PROC0"] = "OUI";

    ListSyntaxMapContainer listeResu;
    SyntaxMapContainer dict2;
    dict2.container["CHAM_GD"] = getName();
    listeResu.push_back( dict2 );
    dict.container["RESU"] = listeResu;

    cmdSt.define( dict );

    ASTERINTEGER op = 39;
    CALL_EXECOP( &op );

    return true;
};

/** @typedef FieldOnCellsReal */
using FieldOnCellsReal = FieldOnCells< ASTERDOUBLE >;
using FieldOnCellsRealPtr = std::shared_ptr< FieldOnCellsReal >;

/** @typedef FieldOnCellsLong */
using FieldOnCellsLong = FieldOnCells< ASTERINTEGER >;
using FieldOnCellsLongPtr = std::shared_ptr< FieldOnCellsLong >;

/** @typedef FieldOnCellsComplex */
using FieldOnCellsComplex = FieldOnCells< ASTERCOMPLEX >;
using FieldOnCellsComplexPtr = std::shared_ptr< FieldOnCellsComplex >;

/** @typedef FieldOnCellsChar8 */
using FieldOnCellsChar8 = FieldOnCells< JeveuxChar8 >;
using FieldOnCellsChar8Ptr = std::shared_ptr< FieldOnCellsChar8 >;

#endif /* FIELDONCELLS_H_ */
