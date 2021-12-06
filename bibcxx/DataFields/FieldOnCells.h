#ifndef FIELDONCELLS_H_
#define FIELDONCELLS_H_

/**
 * @file FieldOnCells.h
 * @brief Fichier entete de la classe FieldOnCells
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

#include <string>
#include <assert.h>

#include "astercxx.h"
#include "aster_fort_superv.h"
#include "aster_fort_ds.h"

#include "Behaviours/BehaviourProperty.h"
#include "Discretization/ElementaryCharacteristics.h"
#include "DataFields/DataField.h"
#include "DataFields/SimpleFieldOnCells.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "PythonBindings/LogicalUnitManager.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/Exceptions.h"


/**
 * @class FieldOnCells
 * @brief Cette classe template permet de definir un champ aux éléments Aster
 * @author Nicolas Sellenet
 */
template < class ValueType > class FieldOnCells : public DataField {
  private:
    typedef SimpleFieldOnCells< ValueType > SimpleFieldOnCellsValueType;
    typedef boost::shared_ptr< SimpleFieldOnCellsValueType > SimpleFieldOnCellsValueTypePtr;

    /** @brief Vecteur Jeveux '.CELD' */
    JeveuxVectorLong _descriptor;
    /** @brief Vecteur Jeveux '.CELK' */
    JeveuxVectorChar24 _reference;
    /** @brief Vecteur Jeveux '.CELV' */
    JeveuxVector< ValueType > _valuesList;
    /** @brief Modele */
    ModelPtr _model;
    /** @brief Finite element description */
    FiniteElementDescriptorPtr _dofDescription;
    /** @brief jeveux vector '.TITR' */
    JeveuxVectorChar80 _title;

  public:
    /**
     * @typedef FieldOnCellsPtr
     * @brief Pointeur intelligent vers un FieldOnCells
     */
    typedef boost::shared_ptr< FieldOnCells > FieldOnCellsPtr;

    /**
     * @brief Constructeur
     * @param name Nom Jeveux du champ aux éléments
     */
    FieldOnCells( const std::string name )
        : DataField( name, "CHAM_ELEM" ),
          _descriptor( JeveuxVectorLong( getName() + ".CELD" ) ),
          _reference( JeveuxVectorChar24( getName() + ".CELK" ) ),
          _valuesList( JeveuxVector< ValueType >( getName() + ".CELV" ) ), _model( nullptr ),
          _title( JeveuxVectorChar80( getName() + ".TITR" ) ){};

    /**
     * @brief Constructeur
     */
    FieldOnCells( )
        : FieldOnCells( DataStructureNaming::getNewName() ){};


    /**
     * @brief Constructeur à partir d'une carte compor et d'un modèle
     * @param model Modèle
     * @param behaviour Carte Compor
     * @param typcham Type de champ à calculer
     */
    FieldOnCells(const ModelPtr &model, const BehaviourPropertyPtr behaviour,
                 const std::string& typcham,const ElementaryCharacteristicsPtr carael = nullptr)
        : FieldOnCells( ){
            std::string inName = getName();
            std::string carele = " ";
            std::string test;
            std::string option;
            std::string nompar;
            test=typcham;
            test.resize(4);

            if (test=="ELGA"){
                option="TOU_INI_ELGA";
            }
            else if(test=="ELNO"){
                option="TOU_INI_ELNO";
            }
            else{
                AS_ASSERT(false)
            };
            if (typcham==test+"_SIEF_R"){
                nompar="PSIEF_R";
            }
            else if(typcham==test+"_VARI_R"){
                nompar="PVARI_R";
            }
            else{
                AS_ASSERT(false)
            };

            if(carael) carele = carael->getName();

            ASTERINTEGER iret = 0;
            auto fed = model->getFiniteElementDescriptor();

            _model = model;
            _dofDescription = fed;
            auto dcel = boost::make_shared<SimpleFieldOnCellsValueType>();
            auto compor = behaviour->getBehaviourField();
            CALLO_CESVAR(carele, compor->getName(), fed->getName(), dcel->getName());
            CALLO_ALCHML(fed->getName(), option, nompar, JeveuxMemoryTypesNames[Permanent],
                         getName(),&iret, dcel->getName());
            AS_ASSERT(iret==0);

            AS_ASSERT(updateValuePointers());
        };


    /**
     * @brief Move constructor
     * @param tocopy FieldOnCells object
     */
    FieldOnCells(FieldOnCells&& toCopy )
        : FieldOnCells(){

          _descriptor = toCopy._descriptor;
          _reference  = toCopy._reference;
          _valuesList = toCopy._valuesList;
          _title      = toCopy._title;
          _dofDescription = toCopy._dofDescription;
          _model = toCopy._model;

          AS_ASSERT(updateValuePointers());
    };

     /**
     * @brief Copy constructor
     */
    FieldOnCells( const std::string &name, const FieldOnCells &toCopy )
        : FieldOnCells( name ){
        // JeveuxVector to be duplicated
        *( _descriptor ) = *( toCopy._descriptor );
        *( _reference ) = *( toCopy._reference );
        *( _valuesList ) = *( toCopy._valuesList );
        *( _title ) = *( toCopy._title );
        // Pointers to be copied
        _dofDescription = toCopy._dofDescription;
        _model = toCopy._model;

        AS_ASSERT(updateValuePointers());
    }


    /**
     * @brief Copy constructor
     */
    FieldOnCells( const FieldOnCells &toCopy )
        : FieldOnCells(DataStructureNaming::getNewName(), toCopy) {};

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
    };

    /**
     * @brief
     * @return
     */
    SimpleFieldOnCellsValueTypePtr exportToSimpleFieldOnCells() {
        SimpleFieldOnCellsValueTypePtr toReturn(
            new SimpleFieldOnCellsValueType() );
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
    ASTERBOOL isSimilarTo(const FieldOnCells<ValueType> &tmp2 ) const {
        bool similar = (_descriptor->size() == tmp2._descriptor->size());
        similar = (similar && (_reference->size() == tmp2._reference->size()));
        similar = (similar && (_valuesList->size() == tmp2._valuesList->size()));
        return similar;
    }

    /**
     * @brief Get the model
     */
    ModelPtr getModel() const {
        if (_model != nullptr && _model->isEmpty() )
            raiseAsterError( "Model is empty" );

        return _model;
    };

    /**
     * @brief Get the mesh
     */
    BaseMeshPtr getMesh() const {
        return _model->getMesh();
    };

    /**
     * @brief Set the description of finite elements
     * @param curDesc object FiniteElementDescriptorPtr
     */
    void setDescription( FiniteElementDescriptorPtr &curDesc ) {
        if ( _dofDescription )
            raiseAsterError( "FiniteElementDescriptor already set" );

        _dofDescription = curDesc;
    };

    /**
     * @brief Definition du modele
     * @param currentMesh objet Model sur lequel la charge reposera
     */
    ASTERBOOL setModel( ModelPtr &currentModel ) {
        if ( currentModel->isEmpty() )
            raiseAsterError( "Model is empty" );

        _model = currentModel;
        return true;
    };

    /**
     * @brief Update field and build FiniteElementDescriptor if necessary
     */
    ASTERBOOL build() {
        if ( _dofDescription == nullptr && updateValuePointers() ) {
            if ( _model == nullptr )
                raiseAsterError( "Model is empty" );

            const std::string ligrel = trim( ( *_reference )[0].toString() );

            if ( ligrel.substr( 0, 8 ) == getName().substr( 0, 8 ) ) {
                _dofDescription =
                    boost::make_shared< FiniteElementDescriptor >( ligrel, getMesh() );
            } else {
                _dofDescription = _model->getFiniteElementDescriptor();
            }
        }
        return true;
    };

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
     */
    ASTERBOOL updateValuePointers() const {
        bool retour = _descriptor->updateValuePointer();
        retour = ( retour && _reference->updateValuePointer() );
        retour = ( retour && _valuesList->updateValuePointer() );
        return retour;
    };


    /**
     * @brief Transormer les valeurs de _valuesList en appliquant
     *         la fonction "func" à chaque valeur
     * @return renvoie un nouveau objet de FieldOnCells
     *         avec les valeurs transformées
     */
    template < class type = ValueType >
    typename std::enable_if<std::is_same< type,ASTERDOUBLE >::value,FieldOnCells<ValueType>>::type
    transform(PyObject* func) {
        if(!PyCallable_Check(func))
            raiseAsterError("Input parameter to the transform \
        method should be a callable Python object");

        FieldOnCells<ValueType> tmp(*this);
        AS_ASSERT( updateValuePointers() );

        ASTERINTEGER size = _valuesList->size();
        for(auto i=0;i<size;i++){
            PyObject* res = PyObject_CallFunction(func, "d", (*_valuesList)[i]);
            if(PyFloat_Check(res)){
                tmp[i] = (ASTERDOUBLE)PyFloat_AsDouble(res);
            }
            else{
                PyErr_Format(PyExc_ValueError, "Returned value of \
                    type different from ASTERDOUBLE");
                PyErr_Print();
            }
            Py_XDECREF(res);
        }

        return tmp;
    };

    template < class type = ValueType >
    typename std::enable_if<std::is_same< type,ASTERCOMPLEX>::value,FieldOnCells<ValueType>>::type
    transform(PyObject* func) {
        if(!PyCallable_Check(func))
            raiseAsterError("Input parameter to the transform \
        method should be a callable Python object");

        FieldOnCells<ValueType> tmp(*this);
        AS_ASSERT( _valuesList->updateValuePointer() );

        ASTERINTEGER size =  _valuesList->size();

        Py_complex val;
        for(auto i=0;i<size;i++){
            val.real = (*_valuesList)[i].real();
            val.imag = (*_valuesList)[i].imag();
            PyObject* res = PyObject_CallFunction(func, "D", val);
            if(PyComplex_Check(res)){
                ASTERDOUBLE re = (ASTERDOUBLE)PyComplex_RealAsDouble(res);
                ASTERDOUBLE im = (ASTERDOUBLE)PyComplex_ImagAsDouble(res);
                tmp[i] = {re,im};
            }else{
                PyErr_Format(PyExc_ValueError, "Returned value of \
                    type different from ASTERCOMPLEX");
                PyErr_Print();
            }
            //Py_DECREF(res);
            Py_XDECREF(res);
        }
        return tmp;
    };


    // OVERLOADING C++ OPERATORS

    /**
     * @brief Unary Minus overloading
     * @return Updated field
     */
    FieldOnCells< ValueType > operator-() const {
        FieldOnCells<ValueType> tmp(*this);
        AS_ASSERT( _valuesList->updateValuePointer() );
        ASTERINTEGER size = _valuesList->size();
        for ( auto pos = 0; pos < size; ++pos ) tmp[pos] = -(*this )[pos];
        return tmp;
    };

    /**
     * @brief Shorthand + operator assignement
     * @return Updated field
     */
    FieldOnCells< ValueType > &operator+=( const FieldOnCells< ValueType > &rhs ) {
        AS_ASSERT( _valuesList->updateValuePointer() );
        AS_ASSERT( rhs.updateValuePointers() );

        if (!this->isSimilarTo(rhs)){
            raiseAsterError("Fields have incompatible shapes");
        }

        ASTERINTEGER size = _valuesList->size();
        for ( int pos = 0; pos < size; ++pos )
            (*this )[pos] += rhs[pos];

        return *this;
    };

    /**
     * @brief Shorthand - operator assignement
     * @return Updated field
     */
    FieldOnCells< ValueType > &operator-=( const FieldOnCells< ValueType > &rhs ) {
        AS_ASSERT( _valuesList->updateValuePointer() );
        AS_ASSERT( rhs.updateValuePointers() );

        if (!this->isSimilarTo(rhs)){
            raiseAsterError("Fields have incompatible shapes");
        }

        ASTERINTEGER size = _valuesList->size();
        for ( int pos = 0; pos < size; ++pos )
            (*this )[pos] -= rhs[pos];

        return *this;
    };

    /**
     * @brief Subscripting overloading
     * @param i subscript
     * @return value at position i
     */
    ValueType &operator[]( int i ) {
#ifdef ASTER_DEBUG_CXX
        AS_ASSERT( 0 <= i && i < this->size() );
#endif
        return _valuesList->operator[]( i );
    };

    const ValueType &operator[]( int i ) const {
      return const_cast<ValueType&>(const_cast<FieldOnCells<ValueType>*>(this)->operator[](i));
    };

    /**
     * @brief Plus overloading
     * @return New field
     */
    FieldOnCells<ValueType> operator+(const FieldOnCells<ValueType>& rhs){
        FieldOnCells<ValueType> tmp(*this);
        AS_ASSERT( _valuesList->updateValuePointer() );
        AS_ASSERT( rhs.updateValuePointers() );

        if (!tmp.isSimilarTo(rhs))
            raiseAsterError("Fields have incompatible shapes");

        ASTERINTEGER size = rhs._valuesList->size();
        for(auto i=0;i<size;i++){
            (*tmp._valuesList)[i] = (*_valuesList)[i] + (*rhs._valuesList)[i];
        }

        return tmp;
    };

    /**
     * @brief Minus overloading
     * @return New field
     */
    FieldOnCells<ValueType> operator-(const FieldOnCells<ValueType>& rhs){
        FieldOnCells<ValueType> tmp(*this);
        AS_ASSERT( _valuesList->updateValuePointer() );
        AS_ASSERT( rhs.updateValuePointers() );
        ASTERINTEGER size = rhs._valuesList->size();

        if (!tmp.isSimilarTo(rhs))
            raiseAsterError("Fields have incompatible shapes");

        for(auto i=0;i<size;i++){
            (*tmp._valuesList)[i] = (*_valuesList)[i] - (*rhs._valuesList)[i];
        }

        return tmp;
    };

    /**
     * @brief Multiply by a scalar on right overloading
     * @return New field
     */

    friend FieldOnCells< ValueType > operator*(const FieldOnCells< ValueType >& lhs,
                                                const ASTERDOUBLE& scal ) {

        AS_ASSERT( lhs.updateValuePointers() )

        ASTERINTEGER taille = lhs._valuesList->size();
        FieldOnCells< ValueType > tmp(lhs);
        for ( int pos = 0; pos < taille; ++pos )
            tmp[pos] = lhs[pos] * scal;

        return tmp;
    };

    /**
     * @brief Multiply by a scalar on left overloading
     * @return New field
     */

    friend FieldOnCells< ValueType > operator*( const ASTERDOUBLE &scal,
                                                FieldOnCells< ValueType >& rhs ) {
        return rhs * scal;
    };


    // some getters

    /**
     * @brief Get values of the field
     *
     */
    const JeveuxVector< ValueType >& getValues( ) const
    {
        return _valuesList;
    }

    /**
     * @brief Set the Values object
     *
     * @param value Value to affect
     */
    void setValues( const ValueType &value ) {
        AS_ASSERT( _valuesList->updateValuePointer() );
        const int taille = _valuesList->size();

        for ( int pos = 0; pos < taille; ++pos )
            ( *this )[pos] = value;
    };

    /**
     * @brief Size of the FieldOnNodes
     */
    const ASTERINTEGER size( ) const
    {
        return _valuesList->size();
    }

    // norm and dot methods


    /**
     * @brief Comput norm
     * @param normType Type of norm ("NORM_1","NORM_2","NORM_INFINITY")
     */
    template < class type = ValueType >
    typename std::enable_if< std::is_same< type, ASTERDOUBLE >::value, ASTERDOUBLE>::type
    norm(const std::string normType) const{
        ASTERDOUBLE norme = 0.0;
        ASTERINTEGER beg = 0, end = 0, nbgrp = 0;

        const int rank = getMPIRank();

        AS_ASSERT( _valuesList->updateValuePointer() );
        AS_ASSERT( _descriptor->updateValuePointer() );

        /*if ( getMesh()->isParallel() ) {
            AS_ASSERT(false);
        }*/

        JeveuxVectorLong CellsRank = getMesh()->getCellsRank();
        bool retour = CellsRank->updateValuePointer();

        if(!_model || _model->isEmpty()){
            raiseAsterError("Model not assigned to the FieldOnCells or empty");
        }

        JeveuxCollectionLong collec = _model->getFiniteElementDescriptor()->getListOfGroupOfCells();
        JeveuxVectorLong descr = _descriptor;
        nbgrp =  (*descr)[1];

        for(auto i = 0; i < nbgrp; i++){

            ASTERINTEGER adress =   (*descr)[4+i];
            if((*descr)[adress + 2] == 0) continue;

            ASTERINTEGER nel = (*descr)[adress];
            auto liel = collec->getObject(i+1);

            if( normType == "NORM_1"){
                for(auto p = 0; p < nel; p++){

                    if((*CellsRank)[liel[p]-1] != rank) continue;
                    beg = (*descr)[adress + 3 + 4 * p + 4] - 1;
                    end = beg + (*descr)[adress + 3 + 4 * p + 3];

                    for( int pos = beg; pos < end; ++pos ){
                        norme += std::abs(( *this )[pos]);
                    }
                }
            }
            else if( normType == "NORM_2"){
                for(auto p = 0; p < nel; p++){

                    if((*CellsRank)[liel[p]-1] != rank) continue;
                    beg = (*descr)[adress + 3 + 4 * p + 4] - 1;
                    end = beg + (*descr)[adress + 3 + 4 * p + 3];

                    for( int pos = beg; pos < end; ++pos ){
                        norme += ( *this )[pos] * ( *this )[pos];
                    }
                }
            }
            else if( normType == "NORM_INFINITY") {
                for(auto p = 0; p < nel; p++){

                  if((*CellsRank)[liel[p]-1] != rank) continue;
                    beg = (*descr)[adress + 3 + 4 * p + 4] - 1;
                    end = beg + (*descr)[adress + 3 + 4 * p + 3];

                    for( int pos = beg; pos < end; ++pos ){
                        norme = std::max(norme, std::abs(( *this )[pos]));
                    }
                }
            }
            else
            {
                AS_ASSERT(false);
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

        if ( normType == "NORM_2" ) norme = std::sqrt( norme );

        return norme;
    }

    /**
     * @brief Dot product
     * @param tmp object FieldOnCellsPtr
     */
    template < class type = ValueType >
    typename std::enable_if< std::is_same< type, ASTERDOUBLE >::value, ASTERDOUBLE>::type
    dot( const FieldOnCellsPtr &tmp ) const {
        AS_ASSERT( tmp->updateValuePointers());
        AS_ASSERT( _valuesList->updateValuePointer() );
        ASTERINTEGER taille = _valuesList->size();

        if( taille != tmp->size() )
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

    try {
        ASTERINTEGER op = 39;
        CALL_EXECOP( &op );
    } catch ( ... ) {
        throw;
    }

    return true;
};

/** @typedef FieldOnCellsReal Class d'une carte de double */
typedef FieldOnCells< ASTERDOUBLE > FieldOnCellsReal;

/**
 * @typedef FieldOnCellsPtrReal
 * @brief Definition d'un champ aux éléments de double
 */
typedef boost::shared_ptr< FieldOnCellsReal > FieldOnCellsRealPtr;

/** @typedef FieldOnCellsLong Class d'une carte de long */
typedef FieldOnCells< ASTERINTEGER > FieldOnCellsLong;

/**
 * @typedef FieldOnCellsPtrLong
 * @brief Definition d'un champ aux éléments de long
 */
typedef boost::shared_ptr< FieldOnCellsLong > FieldOnCellsLongPtr;

/** @typedef FieldOnCellsLong Class d'une carte de complex */
typedef FieldOnCells< ASTERCOMPLEX > FieldOnCellsComplex;

/**
 * @typedef FieldOnCellsPtrComplex
 * @brief Definition d'un champ aux éléments de complexes
 */
typedef boost::shared_ptr< FieldOnCellsComplex > FieldOnCellsComplexPtr;

#endif /* FIELDONCELLS_H_ */
