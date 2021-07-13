#ifndef FIELDONCELLS_H_
#define FIELDONCELLS_H_

/**
 * @file FieldOnCells.h
 * @brief Fichier entete de la classe FieldOnCells
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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
#include "MemoryManager/JeveuxVector.h"
#include "DataFields/DataField.h"
#include "DataFields/SimpleFieldOnCells.h"
#include "Modeling/Model.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Behaviours/BehaviourProperty.h"
#include "PythonBindings/LogicalUnitManager.h"
#include "Supervis/CommandSyntax.h"


/**
 * @class FieldOnCells
 * @brief Cette classe template permet de definir un champ aux éléments Aster
 * @author Nicolas Sellenet
 */
template < class ValueType > class FieldOnCells : public DataField {
  private:
    typedef SimpleFieldOnCells< ValueType > SimpleFieldOnCellsValueType;
    typedef boost::shared_ptr< SimpleFieldOnCellsReal >
        SimpleFieldOnCellsValueTypePtr;

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
     * @param memType Mémoire d'allocation
     */
    FieldOnCells( const JeveuxMemory memType = Permanent )
        : DataField( memType, "CHAM_ELEM" ),
          _descriptor( JeveuxVectorLong( getName() + ".CELD" ) ),
          _reference( JeveuxVectorChar24( getName() + ".CELK" ) ),
          _valuesList( JeveuxVector< ValueType >( getName() + ".CELV" ) ), _model( nullptr ),
          _title( JeveuxVectorChar80( getName() + ".TITR" ) ){};


    /**
     * @brief Constructeur à partir d'une carte compor et d'un modèle
     * @param model Modèle
     * @param behaviour Carte Compor
     * @param typcham Type de champ à calculer
     * @param memType Mémoire d'allocation
     */
    FieldOnCells(const ModelPtr &model, const BehaviourPropertyPtr behaviour,
                 const std::string& typcham, const JeveuxMemory memType = Permanent )
        : FieldOnCells( memType){
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

            ASTERINTEGER iret = 0;
            _model=model;
            _dofDescription=model->getFiniteElementDescriptor();
            auto fed = model->getFiniteElementDescriptor();
            auto dcel = boost::make_shared<SimpleFieldOnCellsValueType>( getMemoryType() );
            auto compor = behaviour->getBehaviourField();
            CALLO_CESVAR(carele, compor->getName(), fed->getName(), dcel->getName());
            CALLO_ALCHML(fed->getName(), option, nompar, JeveuxMemoryTypesNames[getMemoryType()],
                         getName(),&iret, dcel->getName());
            AS_ASSERT(iret==0);};


    /**
     * @brief Move constructor
     * @param tocopy FieldOnCells object
     */
    FieldOnCells(FieldOnCells&& toCopy )
        : DataField( toCopy.getMemoryType(), "CHAM_ELEM" ),
          _descriptor( JeveuxVectorLong( getName() + ".CELD" ) ),
          _reference( JeveuxVectorChar24( getName() + ".CELK" ) ),
          _valuesList( JeveuxVector< ValueType >( getName() + ".CELV" ) ),
          _title( JeveuxVectorChar80( getName() + ".TITR" ) ){

          _descriptor = toCopy._descriptor;
          _reference  = toCopy._reference;
          _valuesList = toCopy._valuesList;
          _title      = toCopy._title;
          _dofDescription = toCopy._dofDescription;
          _model = toCopy._model;
    };

     /**
     * @brief Copy constructor
     */
    FieldOnCells( const std::string &name, const FieldOnCells &toCopy )
        : DataField( name, "CHAM_ELEM", toCopy.getMemoryType() ),
          _descriptor( JeveuxVectorLong( getName() + ".CELD" ) ),
          _reference( JeveuxVectorChar24( getName() + ".CELK" ) ),
          _valuesList( JeveuxVector< ValueType >( getName() + ".CELV" ) ),
          _title( JeveuxVectorChar80( getName() + ".TITR" ) ){
        // JeveuxVector to be duplicated
        *( _descriptor ) = *( toCopy._descriptor );
        *( _reference ) = *( toCopy._reference );
        *( _valuesList ) = *( toCopy._valuesList );
        *( _title ) = *( toCopy._title );
        // Pointers to be copied
        _dofDescription = toCopy._dofDescription;
        _model = toCopy._model;
    }

    /**
     * @brief Wrap of copy constructor
     */
    FieldOnCells duplicate() { return *this; }

    /**
     * @brief Copy constructor
     */
    FieldOnCells( const FieldOnCells &toCopy )
        : FieldOnCells(ResultNaming::getNewResultName(), toCopy) {};

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
            new SimpleFieldOnCellsValueType( getMemoryType() ) );
        const std::string resultName = toReturn->getName();
        const std::string inName = getName();
        const std::string copyNan( "OUI" );
        CALLO_CELCES_WRAP( inName, JeveuxMemoryTypesNames[getMemoryType()], resultName );
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
            throw std::runtime_error( "Model is empty" );
        return _model;
    };

    /**
     * @brief Set the description of finite elements
     * @param curDesc object FiniteElementDescriptorPtr
     */
    void setDescription( FiniteElementDescriptorPtr &curDesc ) {
        if ( _dofDescription )
            throw std::runtime_error( "FiniteElementDescriptor already set" );
        _dofDescription = curDesc;
    };

    /**
     * @brief Definition du modele
     * @param currentMesh objet Model sur lequel la charge reposera
     */
    ASTERBOOL setModel( ModelPtr &currentModel ) {
        if ( currentModel->isEmpty() )
            throw std::runtime_error( "Model is empty" );
        _model = currentModel;
        return true;
    };

    /**
     * @brief Update field and build FiniteElementDescriptor if necessary
     */
    ASTERBOOL build() {
        if ( _dofDescription == nullptr && updateValuePointers() ) {
            if ( _model == nullptr )
                throw std::runtime_error( "Model is empty" );
            _dofDescription = _model->getFiniteElementDescriptor();
        }
        return true;
    };

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
     */
    ASTERBOOL updateValuePointers() {
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
        if(!PyCallable_Check(func)) throw std::runtime_error("Input parameter to the transform \
        method should be a callable Python object");
        FieldOnCells<ValueType> tmp(*this);
        ASTERBOOL ret = true;
        if(_valuesList.isEmpty()) ret = updateValuePointers();
        if(ret){
            ASTERINTEGER size =  _valuesList->size();
            //vect.resize(size);
            for(auto i=0;i<size;i++){
                PyObject* res = PyObject_CallFunction(func, "d", (*_valuesList)[i]);
                if(PyFloat_Check(res)){
                    tmp[i] = (ASTERDOUBLE)PyFloat_AsDouble(res);
                }else{
                    PyErr_Format(PyExc_ValueError, "Returned value of \
                    type different from ASTERDOUBLE");
                    PyErr_Print();
                }
                //Py_DECREF(res);
                Py_XDECREF(res);
            }
            return tmp;
        }
    };

    template < class type = ValueType >
    typename std::enable_if<std::is_same< type,ASTERCOMPLEX>::value,FieldOnCells<ValueType>>::type
    transform(PyObject* func) {
        if(!PyCallable_Check(func)) throw std::runtime_error("Input parameter to the transform \
        method should be a callable Python object");
        FieldOnCells<ValueType> tmp(*this);
        if(_valuesList.isEmpty()) throw std::runtime_error("FieldOnCells of complex type is empty");
        ASTERINTEGER size =  _valuesList->size();
        //vect.resize(size);
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
        //ASTERBOOL ret = updateValuePointers();
        FieldOnCells<ValueType> tmp(*this);
        ASTERINTEGER size = _valuesList->size();
        for ( auto pos = 0; pos < size; ++pos )tmp[pos] = -(*this )[pos];
        return tmp;
    };

    /**
     * @brief Shorthand + operator assignement
     * @return Updated field
     */
    FieldOnCells< ValueType > &operator+=( const FieldOnCells< ValueType > &rhs ) {
        ASTERBOOL ret = true;
        if(_valuesList.isEmpty()) ret = updateValuePointers();
        if (rhs._valuesList.isEmpty()) {
            ret = ret && const_cast<FieldOnCells< ValueType >&> (rhs).updateValuePointers();
        }
        if(ret){
            if (!this->isSimilarTo(rhs)){
                throw std::runtime_error("Fields have incompatible shapes");
            }
            ASTERINTEGER size = _valuesList->size();
            for ( int pos = 0; pos < size; ++pos ) (*this )[pos] = ( *this )[pos] + rhs[pos];
            return *this;
        }else{
            throw  std::runtime_error("Unable to use the operator+= : \
            Maye one of the fieldOnCells objects is empty");
        }
    };

    /**
     * @brief Shorthand - operator assignement
     * @return Updated field
     */
    FieldOnCells< ValueType > &operator-=( const FieldOnCells< ValueType > &rhs ) {
        ASTERBOOL ret = true;
        if(_valuesList.isEmpty()) ret = updateValuePointers();
        if (rhs._valuesList.isEmpty()){
            ret = ret && const_cast<FieldOnCells< ValueType >&> (rhs).updateValuePointers();
        }
        if(ret){
            if (!this->isSimilarTo(rhs)){
                throw std::runtime_error("Fields have incompatible shapes");
            }
            ASTERINTEGER size = _valuesList->size();
            for ( int pos = 0; pos < size; ++pos ) (*this )[pos] = ( *this )[pos] - rhs[pos];
            return *this;
        }else{
            throw  std::runtime_error("Unable to use the operator-= : Maybe \
            one of the fieldOnCells objects is empty");
        }
    };

    /**
     * @brief Subscripting overloading
     * @param i subscript
     * @return value at position i
     */
    ValueType &operator[]( int i ) {
        if( 0 <= i  && i < this->size()){
            return _valuesList->operator[]( i );
        }else{
            throw std::runtime_error("Index out of range");
        }
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
        ASTERBOOL ret = true;
        if(_valuesList.isEmpty()) ret = updateValuePointers();
        if (rhs._valuesList.isEmpty()){
            ret = ret && const_cast<FieldOnCells< ValueType >&> (rhs).updateValuePointers();
        }
        if(ret){
            if (!tmp.isSimilarTo(rhs)) throw std::runtime_error("Fields have incompatible shapes");
            ASTERINTEGER size = rhs._valuesList->size();
             for(auto i=0;i<size;i++){
                (*tmp._valuesList)[i] = (*_valuesList)[i] + (*rhs._valuesList)[i];
             }
             return tmp;
        }else{
            throw  std::runtime_error("Unable to use the operator + : \
            Maybe one of the fieldOnCells objects is empty");
        }
    };

    /**
     * @brief Minus overloading
     * @return New field
     */
    FieldOnCells<ValueType> operator-(const FieldOnCells<ValueType>& rhs){
        FieldOnCells<ValueType> tmp(*this);
        ASTERBOOL ret = true;
        if(_valuesList.isEmpty()) ret = updateValuePointers();
        if (rhs._valuesList.isEmpty()){
            ret = ret && const_cast<FieldOnCells< ValueType >&> (rhs).updateValuePointers();
        }
        ASTERINTEGER size = rhs._valuesList->size();
        if(ret){
            if (!tmp.isSimilarTo(rhs)) throw std::runtime_error("Fields have incompatible shapes");
            for(auto i=0;i<size;i++){
                (*tmp._valuesList)[i] = (*_valuesList)[i] - (*rhs._valuesList)[i];
            }
            return tmp;
        }else{
            throw  std::runtime_error("Unable to use the operator - :  \
            Maybe one of the fieldOnCells objects is empty");
        }
    };

    /**
     * @brief Multiply by a scalar on right overloading
     * @return New field
     */

    friend FieldOnCells< ValueType > operator*(const FieldOnCells< ValueType >& lhs,
                                                const ASTERDOUBLE& scal ) {
        ASTERBOOL ret = true;
        if (lhs._valuesList.isEmpty()){
            ret = const_cast<FieldOnCells< ValueType >&> (lhs).updateValuePointers();
        }
        if(ret){
            ASTERINTEGER taille = lhs._valuesList->size();
            FieldOnCells< ValueType > tmp(lhs);
            for ( int pos = 0; pos < taille; ++pos ) tmp[pos] = lhs[pos] * scal;
            return tmp;
        }else{
            throw  std::runtime_error("Unable to use the operator * :\
             Maybe the fieldOnCells object is empty");
        }
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
        bool retour = _valuesList->updateValuePointer();
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
        ASTERBOOL ret =  _valuesList->updateValuePointer();
        if(ret){
            int taille = _valuesList->size();
            if( normType == "NORM_1"){
                for( int pos = 0; pos < taille; ++pos ){
                    norme += std::abs(( *this )[pos]);
                }
            }
            else if( normType == "NORM_2"){
                for( int pos = 0; pos < taille; ++pos ){
                    norme += ( *this )[pos] * ( *this )[pos];
                }
            }

            else if( normType == "NORM_INFINITY") {
                for( int pos = 0; pos < taille; ++pos ){
                    norme = std::max(norme, std::abs(( *this )[pos]));
                }

        }
        }else{
            throw std::runtime_error("Unable to use norm method: \
             Maybe the FieldOnCells Object is empty");
        }
        // square root for l2 norm
        if ( normType == "NORM_2" )  norme = std::sqrt( norme );

        return norme;
    }


    /**
     * @brief Dot product
     * @param tmp object FieldOnCellsPtr
     */
    template < class type = ValueType >
    typename std::enable_if< std::is_same< type, ASTERDOUBLE >::value, ASTERDOUBLE>::type
    dot( const FieldOnCellsPtr &tmp ) const {
        bool retour = tmp->updateValuePointers();
        retour = ( retour && _valuesList->updateValuePointer() );
        ASTERINTEGER taille = _valuesList->size();

        if ( !retour || taille != tmp->size() )
            throw std::runtime_error( "Incompatible size" );

        ASTERDOUBLE ret = 0.0;
        for ( auto pos = 0; pos < taille; ++pos ) {
                ret += ( *this )[pos] * ( *tmp )[pos];
        }
        return ret;
    }

    bool printMedFile( const std::string fileName ) const;

    friend class FieldBuilder;

};

template < class ValueType >
bool FieldOnCells< ValueType >::printMedFile( const std::string fileName ) const {
    LogicalUnitFile a( fileName, Binary, New );
    int retour = a.getLogicalUnit();
    CommandSyntax cmdSt( "IMPR_RESU" );

    SyntaxMapContainer dict;
    dict.container["FORMAT"] = "MED";
    dict.container["UNITE"] = (ASTERINTEGER)retour;

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
