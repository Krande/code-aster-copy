#ifndef FIELDONNODES_H_
#define FIELDONNODES_H_

/**
 * @file FieldOnNodes.h
 * @brief Header of class for FieldOnNodes
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

#include "aster_fort_superv.h"
#include "astercxx.h"

#include "DataFields/DataField.h"
#include "DataFields/MeshCoordinatesField.h"
#include "DataFields/SimpleFieldOnNodes.h"
#include "MemoryManager/JeveuxAllowedTypes.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/BaseMesh.h"
#include "Numbering/DOFNumbering.h"
#include "Numbering/FieldOnNodesDescription.h"
#include "ParallelUtilities/AsterMPI.h"
#include "PythonBindings/LogicalUnitManager.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/Exceptions.h"
#include "Utilities/Blas.h"

/**
 * @struct AllowedFieldType
 * @brief Structure template permettant de limiter le type instanciable de JeveuxVector
 * @tparam T Type autorise
 */
template < typename T >
struct AllowedFieldType; // undefined for bad types!

template <>
struct AllowedFieldType< ASTERINTEGER > {
    static const unsigned short numTypeJeveux = Integer;
};

template <>
struct AllowedFieldType< ASTERDOUBLE > {
    static const unsigned short numTypeJeveux = Real;
};

template <>
struct AllowedFieldType< ASTERCOMPLEX > {
    static const unsigned short numTypeJeveux = Complex;
};

template <>
struct AllowedFieldType< JeveuxChar8 > {
    static const unsigned short numTypeJeveux = Char8;
};

class FieldBuilder;

/**
 * @class FieldOnNodes
 * @brief Template class for FieldOnnodes
 */
template < class ValueType >
class FieldOnNodes : public DataField, private AllowedFieldType< ValueType > {
  private:
    typedef SimpleFieldOnNodes< ValueType > SimpleFieldOnNodesValueType;
    typedef boost::shared_ptr< SimpleFieldOnNodesValueType > SimpleFieldOnNodesValueTypePtr;

    /** @brief Vecteur Jeveux '.DESC' */
    JeveuxVectorLong _descriptor;
    /** @brief Vecteur Jeveux '.REFE' */
    JeveuxVectorChar24 _reference;
    /** @brief Vecteur Jeveux '.VALE' */
    JeveuxVector< ValueType > _valuesList;
    /** @brief Dof description */
    FieldOnNodesDescriptionPtr _dofDescription;
    /** @brief Support mesh */
    BaseMeshPtr _mesh;

  public:
    /** @typedef FieldOnNodesPtr */
    typedef boost::shared_ptr< FieldOnNodes > FieldOnNodesPtr;

    /**
     * @brief Constructor
     * @param name Jeveux name of the field
     */
    FieldOnNodes( const std::string name )
        : DataField( name, "CHAM_NO" ),
          _descriptor( JeveuxVectorLong( getName() + ".DESC" ) ),
          _reference( JeveuxVectorChar24( getName() + ".REFE" ) ),
          _valuesList( JeveuxVector< ValueType >( getName() + ".VALE" ) ),
          _dofDescription( nullptr ),
          _mesh( nullptr ){};

    /** @brief Constructor with automatic name */
    FieldOnNodes() : FieldOnNodes( DataStructureNaming::getNewName() ){};

    /** @brief Copy constructor */
    FieldOnNodes( const std::string &name, const FieldOnNodes &toCopy ) : FieldOnNodes( name ) {
        // JeveuxVector to be duplicated
        *( _descriptor ) = *( toCopy._descriptor );
        *( _reference ) = *( toCopy._reference );
        *( _valuesList ) = *( toCopy._valuesList );
        *( _title ) = *( toCopy._title );
        // Pointers to be copied
        _dofDescription = toCopy._dofDescription;
        _mesh = toCopy._mesh;
    }

    /** @brief Move constructor */
    FieldOnNodes( FieldOnNodes &&other ) : DataField{ std::move( other ) } {
        // Pointers to be moved
        _descriptor = other._descriptor;
        _reference = other._reference;
        _valuesList = other._valuesList;
        _title = other._title;
        _dofDescription = other._dofDescription;
        _mesh = other._mesh;
    }

    /**
     * @brief Copy constructor
     */
    FieldOnNodes( const FieldOnNodes &toCopy )
        : FieldOnNodes( DataStructureNaming::getNewName(), toCopy ){};

    /**
     * @brief Constructor with DOFNumbering
     */
    FieldOnNodes( const BaseDOFNumberingPtr &dofNum ) : FieldOnNodes() {

        _dofDescription = dofNum->getDescription();
        _mesh = dofNum->getMesh();

        const auto intType = AllowedFieldType< ValueType >::numTypeJeveux;
        CALLO_VTCREB_WRAP( getName(), JeveuxMemoryTypesNames[Permanent], JeveuxTypesNames[intType],
                           dofNum->getName() );
    };

    /**
     * @brief Wrap of copy constructor
     * @return new field, copy of the calling field
     */
    FieldOnNodes duplicate() { return *this; }

    /**
     * @brief Constructeur from a MeshCoordinatesFieldPtr&
     */
    FieldOnNodes( MeshCoordinatesFieldPtr &toCopy )
        : DataField( "CHAM_NO" ),
          _descriptor( toCopy->getDescriptor() ),
          _reference( toCopy->getReference() ),
          _valuesList( toCopy->getValues() ),
          _dofDescription( nullptr ),
          _mesh( nullptr ){};

    /**
     * @brief Surcharge de l'operateur []
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    ValueType &operator[]( ASTERINTEGER i ) { return _valuesList->operator[]( i ); };

    /**
     * @brief Surcharge de l'operateur []
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    const ValueType &operator[]( ASTERINTEGER i ) const { return _valuesList->operator[]( i ); };

    /**
     * @brief Check if fields are OK for +, +=, ...
     * @return true if compatible
     */
    bool isSimilarTo( const FieldOnNodes< ValueType > &tmp2 ) {
        CALL_JEMARQ();
        bool similar = ( ( *_descriptor ) == ( *tmp2._descriptor ) );
        similar = ( similar && ( this->_reference->size() == tmp2._reference->size() ) );
        similar = ( similar && ( this->_valuesList->size() == tmp2._valuesList->size() ) );
        similar = ( similar && ( _mesh == tmp2._mesh ) );
        if ( similar ) {
            _descriptor->updateValuePointer();
            tmp2._descriptor->updateValuePointer();
            if ( ( *_descriptor )[1] > 0 || ( *tmp2._descriptor )[1] > 0 )
                similar = ( similar && ( ( *_dofDescription ) == ( *tmp2._dofDescription ) ) );
        }
        CALL_JEDEMA();
        return similar;
    }

    /**
     * @brief PlusEqual overloading
     * @return Updated field
     */
    FieldOnNodes< ValueType > &operator+=( FieldOnNodes< ValueType > const &rhs ) {
        if ( !this->isSimilarTo( rhs ) )
            raiseAsterError( "Fields have incompatible shapes" );
        ( *_valuesList ) += ( *rhs._valuesList );
        return *this;
    };

    /**
     * @brief MinusEqual overloading
     * @return Updated field
     * @todo ajouter une vérification sur la structure des champs
     */
    FieldOnNodes< ValueType > &operator-=( FieldOnNodes< ValueType > const &rhs ) {
        if ( !this->isSimilarTo( rhs ) )
            raiseAsterError( "Fields have incompatible shapes" );
        ( *_valuesList ) -= ( *rhs._valuesList );
        return *this;
    };

    /**
     * @brief TimesEqual overloading
     * @return Updated field
     */
    FieldOnNodes< ValueType > &operator*=( const ASTERDOUBLE &scal ) {
        // AsterBLAS::scal(scal, _valuesList);

        ( *_valuesList ) *= scal;

        return *this;
    };

    /**
     * @brief Unary Minus overloading
     * @return Updated field
     */
    FieldOnNodes< ValueType > operator-() {
        FieldOnNodes< ValueType > tmp( *this );
        ( *tmp._valuesList ) *= ValueType( -1 );
        return tmp;
    };

    // NOTE ON FOLLOWING OPERATOR OVERLAODING
    // friends methods defined inside class body are inline and are hidden from non-ADL lookup
    // passing lhs by value helps optimize chained a+b+c
    // otherwise, both parameters may be const references
    // return the result by value in order to use copy constructor

    /**
     * @brief Multiply by a scalar on right overloading
     * @return New field
     */
    friend FieldOnNodes< ValueType > operator*( FieldOnNodes< ValueType > lhs,
                                                const ASTERDOUBLE &scal ) {

        ( *lhs._valuesList ) *= scal;
        return lhs;
    };

    /**
     * @brief Multiply by a scalar on left overloading
     * @return New field
     */
    friend FieldOnNodes< ValueType > operator*( const ASTERDOUBLE &scal,
                                                FieldOnNodes< ValueType > &rhs ) {
        return rhs * scal;
    };

    /**
     * @brief Plus overloading
     * @return New field
     */
    friend FieldOnNodes< ValueType > operator+( FieldOnNodes< ValueType > lhs,
                                                const FieldOnNodes< ValueType > &rhs ) {
        if ( !lhs.isSimilarTo( rhs ) )
            raiseAsterError( "Fields have incompatible shapes" );
        lhs += rhs;
        return lhs;
    };

    /**
     * @brief Minus overloading
     * @return New field
     * @todo ajouter une vérification sur la structure des champs
     */
    friend FieldOnNodes< ValueType > operator-( FieldOnNodes< ValueType > lhs,
                                                const FieldOnNodes< ValueType > &rhs ) {
        if ( !lhs.isSimilarTo( rhs ) )
            raiseAsterError( "Fields have incompatible shapes" );
        lhs -= rhs;
        return lhs;
    };

    /**
     * @brief Allouer un champ au noeud à partir d'un autre
     * @return renvoit true
     */
    bool allocateFrom( const FieldOnNodes< ValueType > &tmp ) {
        this->_descriptor->deallocate();
        this->_reference->deallocate();
        this->_valuesList->deallocate();

        this->_descriptor->allocate( tmp._descriptor->size() );
        this->_reference->allocate( tmp._reference->size() );
        this->_valuesList->allocate( tmp._valuesList->size() );
        return true;
    };

    /**
     * @brief Renvoit un champ aux noeuds simple (carré de taille nb_no*nbcmp)
     * @return SimpleFieldOnNodesValueTypePtr issu du FieldOnNodes
     */
    SimpleFieldOnNodesValueTypePtr exportToSimpleFieldOnNodes() {
        SimpleFieldOnNodesValueTypePtr toReturn( new SimpleFieldOnNodesValueType() );
        const std::string resultName = toReturn->getName();
        const std::string inName = getName();
        CALLO_CNOCNS_WRAP( inName, JeveuxMemoryTypesNames[Permanent], resultName );
        toReturn->updateValuePointers();
        return toReturn;
    };

    /**
     * @brief Get mesh
     */
    BaseMeshPtr getMesh() const { return _mesh; };

    bool printMedFile( const std::string fileName, bool local = true ) const;

    /**
     * @brief Set the Values object
     *
     * @param value Value to affect
     */
    void setValues( const ValueType &value ) {
        CALL_JEMARQ();
        _valuesList->updateValuePointer();
        _valuesList->assign( value );
        CALL_JEDEMA();
    };

    /**
     * @brief Get values of the field
     *
     */
    const JeveuxVector< ValueType > &getValues() const { return _valuesList; }

    /**
     * @brief Set FieldOnNodes description
     * @param desc object FieldOnNodesDescriptionPtr
     */
    void setDescription( const FieldOnNodesDescriptionPtr &desc ) {
        if ( !desc )
            raiseAsterError( "Empty FieldOnNodesDescription" );
        if ( _dofDescription && _dofDescription->getName() != desc->getName() )
            raiseAsterError( "FieldOnNodesDescription inconsistents" );
        _dofDescription = desc;
    };

    /**
     * @brief Set mesh
     * @param mesh object BaseMeshPtr
     */
    void setMesh( const BaseMeshPtr &mesh ) {
        if ( !mesh )
            raiseAsterError( "Empty Mesh" );

        if ( _mesh && mesh->getName() != _mesh->getName() )
            raiseAsterError( "Meshes inconsistents" );
        _mesh = mesh;
    };

    /**
     * @brief Comput norm
     * @param normType Type of norm ("NORM_1","NORM_2","NORM_INFINITY")
     */
    template < class type = ValueType >
    typename std::enable_if< std::is_same< type, ASTERDOUBLE >::value, ASTERDOUBLE >::type
    norm( const std::string normType ) const {
        CALL_JEMARQ();
        ASTERDOUBLE norme = 0.0;
        _valuesList->updateValuePointer();
        auto taille = _valuesList->size();
        const auto rank = getMPIRank();
        if ( !_mesh )
            raiseAsterError( "Mesh is empty" );
        const JeveuxVectorLong nodesRank = _mesh->getNodesRank();
        nodesRank->updateValuePointer();

        if ( !_dofDescription )
            raiseAsterError( "Description is empty" );
        const VectorLong nodesId = _dofDescription->getNodesFromDOF();

        if ( normType == "NORM_1" ) {
            for ( auto pos = 0; pos < taille; ++pos ) {
                const auto node_id = std::abs( nodesId[pos] );
                if ( ( *nodesRank )[node_id - 1] == rank )
                    norme += std::abs( ( *this )[pos] );
            }
        } else if ( normType == "NORM_2" ) {
            for ( auto pos = 0; pos < taille; ++pos ) {
                const auto node_id = std::abs( nodesId[pos] );
                if ( ( *nodesRank )[node_id - 1] == rank )
                    norme += ( *this )[pos] * ( *this )[pos];
            }
        } else if ( normType == "NORM_INFINITY" ) {
            for ( auto pos = 0; pos < taille; ++pos ) {
                const auto node_id = std::abs( nodesId[pos] );
                if ( ( *nodesRank )[node_id - 1] == rank )
                    norme = std::max( norme, std::abs( ( *this )[pos] ) );
            }
        } else
            raiseAsterError( "Unknown norm" );

#ifdef ASTER_HAVE_MPI
        if ( _mesh->isParallel() ) {
            ASTERDOUBLE norm2 = norme;
            if ( normType == "NORM_1" || normType == "NORM_2" )
                AsterMPI::all_reduce( norm2, norme, MPI_SUM );
            else
                AsterMPI::all_reduce( norm2, norme, MPI_MAX );
        }
#endif

        if ( normType == "NORM_2" )
            norme = std::sqrt( norme );
        CALL_JEDEMA();
        return norme;
    };

    /**
     * @brief Dot product
     * @param tmp object FieldOnNodesDescriptionPtr
     */
    template < class type = ValueType >
    typename std::enable_if< std::is_same< type, ASTERDOUBLE >::value, ASTERDOUBLE >::type
    dot( const FieldOnNodesPtr &tmp ) const {
        CALL_JEMARQ();
        tmp->updateValuePointers();
        _valuesList->updateValuePointer();
        const auto taille = _valuesList->size();

        if ( taille != tmp->size() )
            raiseAsterError( "Incompatible size" );

        if ( !_mesh )
            raiseAsterError( "Mesh is empty" );
        JeveuxVectorLong nodesRank = _mesh->getNodesRank();
        nodesRank->updateValuePointer();

        if ( !_dofDescription )
            raiseAsterError( "Description is empty" );
        const VectorLong nodesId = _dofDescription->getNodesFromDOF();

        const auto rank = getMPIRank();

        ASTERDOUBLE ret = 0.0;
        for ( auto pos = 0; pos < taille; ++pos ) {
            const auto node_id = std::abs( nodesId[pos] );
            if ( ( *nodesRank )[node_id - 1] == rank )
                ret += ( *this )[pos] * ( *tmp )[pos];
        }

#ifdef ASTER_HAVE_MPI
        if ( _mesh->isParallel() ) {
            ASTERDOUBLE ret2 = ret;
            AsterMPI::all_reduce( ret2, ret, MPI_SUM );
        }
#endif
        CALL_JEDEMA();
        return ret;
    }

    bool hasConstantProfile() const {
        CALL_JEMARQ();
        _descriptor->updateValuePointer();
        bool ret = ( ( *_descriptor )[1] <= 0 );
        CALL_JEDEMA();
        return ret;
    }

    /**
     * @brief Size of the FieldOnNodes
     */
    ASTERINTEGER size( void ) const { return _valuesList->size(); }

    /**
     * @brief Get FieldOnNodesDescription
     */
    FieldOnNodesDescriptionPtr getDescription( void ) { return _dofDescription; };

    /**
     * @brief Update field and build FieldOnNodesDescription if necessary
     */
    bool build() {
        if ( !_dofDescription ) {
            CALL_JEMARQ();

            typedef FieldOnNodesDescription FONDesc;
            typedef FieldOnNodesDescriptionPtr FONDescP;

            _reference->updateValuePointer();
            const std::string name2 = trim( ( *_reference )[1].toString() );
            _dofDescription = boost::make_shared< FONDesc >( name2 );
            CALL_JEDEMA();
        }
        return true;
    };

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
     */
    void updateValuePointers() {
        _descriptor->updateValuePointer();
        _reference->updateValuePointer();
        _valuesList->updateValuePointer();
    };

    friend class FieldBuilder;
};

template < class ValueType >
bool FieldOnNodes< ValueType >::printMedFile( const std::string fileName, bool local ) const {
    LogicalUnitFile a( fileName, Binary, New );
    auto retour = a.getLogicalUnit();
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

/** @typedef FieldOnNodesReal */
using FieldOnNodesReal = FieldOnNodes< ASTERDOUBLE >;
using FieldOnNodesRealPtr = boost::shared_ptr< FieldOnNodesReal >;

/** @typedef FieldOnNodesLong */
using FieldOnNodesLong = FieldOnNodes< ASTERINTEGER >;
using FieldOnNodesLongPtr = boost::shared_ptr< FieldOnNodesLong >;

/** @typedef FieldOnNodesComplex*/
using FieldOnNodesComplex = FieldOnNodes< ASTERCOMPLEX >;
using FieldOnNodesComplexPtr = boost::shared_ptr< FieldOnNodesComplex >;

/** @typedef FieldOnNodesChar8 */
using FieldOnNodesChar8 = FieldOnNodes< JeveuxChar8 >;
using FieldOnNodesChar8Ptr = boost::shared_ptr< FieldOnNodesChar8 >;

#endif /* FIELDONNODES_H_ */
