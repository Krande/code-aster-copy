#ifndef FIELDONNODES_H_
#define FIELDONNODES_H_

/**
 * @file FieldOnNodes.h
 * @brief Fichier entete de la classe FieldOnNodes
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

#include <assert.h>

#include "astercxx.h"
#include "aster_fort_superv.h"

#include "DataFields/DataField.h"
#include "DataFields/MeshCoordinatesField.h"
#include "DataFields/SimpleFieldOnNodes.h"
#include "MemoryManager/JeveuxAllowedTypes.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/BaseMesh.h"
#include "Numbering/DOFNumbering.h"
#include "Numbering/FieldOnNodesDescription.h"
#include "ParallelUtilities/MPIContainerUtilities.h"
#include "ParallelUtilities/MPIInfos.h"
#include "PythonBindings/LogicalUnitManager.h"


/**
 * @struct AllowedFieldType
 * @brief Structure template permettant de limiter le type instanciable de JeveuxVectorClass
 * @tparam T Type autorise
 */
template < typename T > struct AllowedFieldType; // undefined for bad types!

template <> struct AllowedFieldType< ASTERINTEGER > {
    static const unsigned short numTypeJeveux = Integer;
};

template <> struct AllowedFieldType< double > { static const unsigned short numTypeJeveux = Real; };

template <> struct AllowedFieldType< RealComplex > {
    static const unsigned short numTypeJeveux = Complex;
};

class FieldBuilder;

/**
 * @class FieldOnNodesClass
 * @brief Cette classe template permet de definir un champ aux noeuds Aster
 * @author Nicolas Sellenet
 */
template < class ValueType >
class FieldOnNodesClass : public DataFieldClass, private AllowedFieldType< ValueType > {
  private:
    typedef SimpleFieldOnNodesClass< ValueType > SimpleFieldOnNodesValueTypeClass;
    typedef boost::shared_ptr< SimpleFieldOnNodesValueTypeClass > SimpleFieldOnNodesValueTypePtr;

    /** @brief Vecteur Jeveux '.DESC' */
    JeveuxVectorLong _descriptor;
    /** @brief Vecteur Jeveux '.REFE' */
    JeveuxVectorChar24 _reference;
    /** @brief Vecteur Jeveux '.VALE' */
    JeveuxVector< ValueType > _valuesList;
    /** @brief Numbering of a FieldOnNodes */
    BaseDOFNumberingPtr _dofNum;
    /** @brief Dof description */
    FieldOnNodesDescriptionPtr _dofDescription;
    /** @brief Support mesh */
    BaseMeshPtr _mesh;
    /** @brief jeveux vector '.TITR' */
    JeveuxVectorChar80 _title;

  public:
    /**
     * @typedef FieldOnNodesPtr
     * @brief Smart pointer to a FieldOnNodes
     */
    typedef boost::shared_ptr< FieldOnNodesClass > FieldOnNodesPtr;

    /**
     * @brief Constructor
     * @param name Jeveux name of the field on nodes
     */
    FieldOnNodesClass( const std::string name )
        : DataFieldClass( name, "CHAM_NO" ), _descriptor( JeveuxVectorLong( getName() + ".DESC" ) ),
          _reference( JeveuxVectorChar24( getName() + ".REFE" ) ),
          _valuesList( JeveuxVector< ValueType >( getName() + ".VALE" ) ), _dofNum( nullptr ),
          _dofDescription( nullptr ), _title( JeveuxVectorChar80( getName() + ".TITR" ) ),
          _mesh( nullptr ){};

    /**
     * @brief Constructor
     * @param memType Type of memory allocation
     */
    FieldOnNodesClass( const JeveuxMemory memType = Permanent )
        : DataFieldClass( memType, "CHAM_NO" ),
          _descriptor( JeveuxVectorLong( getName() + ".DESC" ) ),
          _reference( JeveuxVectorChar24( getName() + ".REFE" ) ),
          _valuesList( JeveuxVector< ValueType >( getName() + ".VALE" ) ), _dofNum( nullptr ),
          _mesh( nullptr ), _dofDescription( nullptr ),
          _title( JeveuxVectorChar80( getName() + ".TITR" ) ){};

    /**
     * @brief Copy constructor
     */
    FieldOnNodesClass( const FieldOnNodesClass &toCopy )
        :DataFieldClass( toCopy.getMemoryType(), "CHAM_NO" ),
          _descriptor( JeveuxVectorLong( getName() + ".DESC" ) ),
          _reference( JeveuxVectorChar24( getName() + ".REFE" ) ),
          _valuesList( JeveuxVector< ValueType >( getName() + ".VALE" ) ), _dofNum( nullptr ),
          _mesh( nullptr ), _dofDescription( nullptr ),
          _title( JeveuxVectorChar80( getName() + ".TITR" ) ) {
        // JeveuxVector to be duplicated
        *(_descriptor) = *(toCopy._descriptor);
        *(_reference) = *(toCopy._reference);
        *(_valuesList) = *(toCopy._valuesList);
        *(_title) = *(toCopy._title);
        // Pointers to be copied
        _dofNum = toCopy._dofNum;
        _dofDescription = toCopy._dofDescription;
        _mesh = toCopy._mesh;
    }

    /**
     * @brief Constructor with DOFNumbering
     */
    FieldOnNodesClass( const BaseDOFNumberingPtr &dofNum, JeveuxMemory memType = Permanent )
        : DataFieldClass( memType, "CHAM_NO" ),
          _descriptor( JeveuxVectorLong( getName() + ".DESC" ) ),
          _reference( JeveuxVectorChar24( getName() + ".REFE" ) ),
          _valuesList( JeveuxVector< ValueType >( getName() + ".VALE" ) ), _dofNum( dofNum ),
          _dofDescription( dofNum->getDescription() ), _mesh( dofNum->getMesh() ),
          _title( JeveuxVectorChar80( getName() + ".TITR" ) )
          {
                if ( !_dofNum )
                    throw std::runtime_error( "DOFNumering is empty" );
                const int intType = AllowedFieldType< ValueType >::numTypeJeveux;
                CALLO_VTCREB_WRAP( getName(), JeveuxMemoryTypesNames[getMemoryType()],
                                JeveuxTypesNames[intType], _dofNum->getName() );
    };

    /**
     * @brief Wrap of copy constructor
     */
    FieldOnNodesClass duplicate() {
        return *this;
    }

    /**
     * @brief Constructeur from a MeshCoordinatesFieldPtr&
     */
    FieldOnNodesClass( MeshCoordinatesFieldPtr &toCopy )
        : DataFieldClass( toCopy->getMemoryType(), "CHAM_NO" ), _descriptor( toCopy->_descriptor ),
          _reference( toCopy->_reference ), _valuesList( toCopy->_valuesList ), _dofNum( nullptr ),
          _dofDescription( nullptr ), _title( JeveuxVectorChar80( getName() + ".TITR" ) ),
          _mesh( nullptr ){};

    /**
     * @brief Surcharge de l'operateur []
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    ValueType &operator[]( int i ) { return _valuesList->operator[]( i ); };

    /**
     * @brief Surcharge de l'operateur []
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    const ValueType &operator[]( int i ) const { return _valuesList->operator[]( i ); };

    /**
     * @brief Check if fields are OK for +, +=, ...
     * @return true if compatible
     */
    bool isSimilarTo(const FieldOnNodesClass< ValueType >  &tmp2 ) const {
        bool similar = (this->_descriptor->size() == tmp2._descriptor->size());
        similar = (similar && (this->_reference->size() == tmp2._reference->size()));
        similar = (similar && (this->_valuesList->size() == tmp2._valuesList->size()));
        return similar;
    }

    /**
     * @brief PlusEqual overloading
     * @return Updated field
     */
    FieldOnNodesClass< ValueType > &operator+=( FieldOnNodesClass< ValueType > const &rhs ) {
        if (!this->isSimilarTo(rhs)) throw std::runtime_error("Fields have incompatible shapes");
        const_cast<FieldOnNodesClass< ValueType >&> (rhs).updateValuePointers() ;
        bool retour = _valuesList->updateValuePointer();
        int taille = _valuesList->size();
        for ( int pos = 0; pos < taille; ++pos )
            ( *this )[pos] = ( *this )[pos] + rhs[pos];
        return *this;
    };

    /**
     * @brief MinusEqual overloading
     * @return Updated field
     * @todo ajouter une vérification sur la structure des champs
     */
    FieldOnNodesClass< ValueType > &operator-=( FieldOnNodesClass< ValueType > const &rhs ) {
        if (!this->isSimilarTo(rhs)) throw std::runtime_error("Fields have incompatible shapes");
        const_cast<FieldOnNodesClass< ValueType >&> (rhs).updateValuePointers() ;
        bool retour = _valuesList->updateValuePointer();
        int taille = _valuesList->size();
        for ( int pos = 0; pos < taille; ++pos )
            ( *this )[pos] = ( *this )[pos] - rhs[pos];
        return *this;
    };

    /**
     * @brief TimesEqual overloading
     * @return Updated field
     */
    FieldOnNodesClass< ValueType > &operator*=( const double &scal ) {
        bool retour = _valuesList->updateValuePointer();
        int taille = _valuesList->size();
        for ( int pos = 0; pos < taille; ++pos )
            ( *this )[pos] = ( *this )[pos] * scal;
        return *this;
    };

    /**
     * @brief Unary Minus overloading
     * @return Updated field
     */
    FieldOnNodesClass< ValueType > &operator-() {
        bool retour = _valuesList->updateValuePointer();
        int taille = _valuesList->size();
        for ( int pos = 0; pos < taille; ++pos )
            ( *this )[pos] = -( *this )[pos];
        return *this;
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
    friend FieldOnNodesClass< ValueType > operator*( FieldOnNodesClass< ValueType > lhs,
                                                     const double &scal ) {
        bool retour = lhs.updateValuePointers();
        int taille = lhs._valuesList->size();
        for ( int pos = 0; pos < taille; ++pos )
            lhs[pos] = lhs[pos] * scal;
        return lhs;
    };

    /**
     * @brief Multiply by a scalar on left overloading
     * @return New field
     */
    friend FieldOnNodesClass< ValueType > operator*( const double &scal,
                                                     FieldOnNodesClass< ValueType > rhs) {
        return rhs * scal;
    };

    /**
     * @brief Plus overloading
     * @return New field
     */
    friend FieldOnNodesClass< ValueType > operator+(FieldOnNodesClass< ValueType > lhs,
                                                    const FieldOnNodesClass< ValueType > &rhs ) {
        if (!lhs.isSimilarTo(rhs)) throw std::runtime_error("Fields have incompatible shapes");
        lhs += rhs;
        return lhs;
    };

    /**
     * @brief Minus overloading
     * @return New field
     * @todo ajouter une vérification sur la structure des champs
     */
    friend FieldOnNodesClass< ValueType > operator-(FieldOnNodesClass< ValueType > lhs,
                                                    const FieldOnNodesClass< ValueType > &rhs ) {
        if (!lhs.isSimilarTo(rhs)) throw std::runtime_error("Fields have incompatible shapes");
        lhs -= rhs;
        return lhs;
    };

    /**
     * @brief Allouer un champ au noeud à partir d'un autre
     * @return renvoit true
     */
    bool allocateFrom( const FieldOnNodesClass< ValueType > &tmp ) {
        this->_descriptor->deallocate();
        this->_reference->deallocate();
        this->_valuesList->deallocate();

        this->_descriptor->allocate( getMemoryType(), tmp._descriptor->size() );
        this->_reference->allocate( getMemoryType(), tmp._reference->size() );
        this->_valuesList->allocate( getMemoryType(), tmp._valuesList->size() );
        return true;
    };

    /**
     * @brief Renvoit un champ aux noeuds simple (carré de taille nb_no*nbcmp)
     * @return SimpleFieldOnNodesValueTypePtr issu du FieldOnNodes
     */
    SimpleFieldOnNodesValueTypePtr exportToSimpleFieldOnNodes() {
        SimpleFieldOnNodesValueTypePtr toReturn(
            new SimpleFieldOnNodesValueTypeClass( getMemoryType() ) );
        const std::string resultName = toReturn->getName();
        const std::string inName = getName();
        CALLO_CNOCNS( inName, JeveuxMemoryTypesNames[getMemoryType()], resultName );
        toReturn->updateValuePointers();
        return toReturn;
    };

    /**
     * @brief Get mesh
     */
    BaseMeshPtr getMesh() const { return _mesh; };

    bool printMedFile( const std::string fileName ) const;

    /**
     * @brief Set DOFNumering
     */
    void setDOFNumbering( const BaseDOFNumberingPtr &dofNum ) {
        if ( _dofNum )
            throw std::runtime_error( "DOFNumbering already set" );
        _dofNum = dofNum;
        _dofDescription = dofNum->getDescription();
        if ( _mesh != nullptr ) {
            const auto name1 = _mesh->getName();
            const auto name2 = _dofNum->getMesh()->getName();
            if ( name1 != name2 )
                throw std::runtime_error( "Meshes inconsistents" );
        }
        else {
            _mesh = dofNum->getMesh();
        }
    };

    /**
     * @brief Set the Values object
     *
     * @param value Value to affect
     */
    void setValues(const ValueType& value)
    {
         bool retour = _valuesList->updateValuePointer();
        const int taille = _valuesList->size();

        for( int pos = 0; pos < taille; ++pos )
            ( *this )[pos] = value;
    };

    /**
     * @brief Get values of the field
     *
     */
    const JeveuxVector< ValueType >& getValues( ) const
    {
        return _valuesList;
    }

    /**
     * @brief Set FieldOnNodes description
     * @param desc object FieldOnNodesDescriptionPtr
     */
    void setDescription( const FieldOnNodesDescriptionPtr &desc ) {
        if ( _dofDescription )
            throw std::runtime_error( "FieldOnNodesDescription already set" );
        _dofDescription = desc;
    };

    /**
     * @brief Set mesh
     * @param mesh object BaseMeshPtr
     */
    void setMesh( const BaseMeshPtr &mesh ) {
        if ( _mesh )
            throw std::runtime_error( "Mesh already set" );
        _mesh = mesh;
        if ( _dofNum != nullptr ) {
            const auto name1 = _mesh->getName();
            const auto name2 = _dofNum->getMesh()->getName();
            if ( name1 != name2 )
                throw std::runtime_error( "Meshes inconsistents" );
        }
    };

    /**
     * @brief Comput norm
     * @param normType Type of norm ("NORM_1","NORM_2","NORM_INFINITY")
     */
    double norm(const std::string normType) const
    {
        double norme = 0.0;
        bool retour =  _valuesList->updateValuePointer();
        int taille = _valuesList->size();
        const int rank = getMPIRank();
        if ( !_mesh )
            throw std::runtime_error( "Mesh is empty" );
        const JeveuxVectorLong nodesRank = _mesh->getNodesRank();
        retour = nodesRank->updateValuePointer();

        if ( !_dofDescription )
            throw std::runtime_error( "Description is empty" );
        const VectorLong nodesId = _dofDescription->getNodesFromDOF();

        if( normType == "NORM_1")
        {
            for( int pos = 0; pos < taille; ++pos )
            {
                const int node_id = std::abs(nodesId[pos]);
                if( (*nodesRank)[node_id-1] == rank )
                    norme += std::abs(( *this )[pos]);
            }
        }
        else if( normType == "NORM_2")
        {
            for( int pos = 0; pos < taille; ++pos )
            {
                const int node_id = std::abs(nodesId[pos]);
                if( (*nodesRank)[node_id-1] == rank )
                    norme += ( *this )[pos] * ( *this )[pos];
            }
        }
        else if( normType == "NORM_INFINITY")
        {
            for( int pos = 0; pos < taille; ++pos )
            {
                const int node_id = std::abs(nodesId[pos]);
                if( (*nodesRank)[node_id-1] == rank )
                    norme = std::max(norme, std::abs(( *this )[pos]));
            }
        }
        else
            throw std::runtime_error( "Unknown norm" );

#ifdef ASTER_HAVE_MPI
        if( _mesh->isParallel() )
        {
            double norm2 = norme;
            MPIContainerUtilities utils;
            if( normType == "NORM_1" || normType == "NORM_2")
                utils.all_reduce(norm2, norme, MPI_SUM);
            else
                utils.all_reduce(norm2, norme, MPI_MAX);
        }
#endif

        if( normType == "NORM_2")
            norme = std::sqrt(norme);

        return norme;
    };

    /**
     * @brief Dot product
     * @param tmp object FieldOnNodesDescriptionPtr
     */
    double dot( const FieldOnNodesPtr &tmp ) const
    {
        bool retour = tmp->updateValuePointers();
        retour = ( retour && _valuesList->updateValuePointer() );
        const int taille = _valuesList->size();

        if( !retour || taille != tmp->size())
            throw std::runtime_error( "Incompatible size" );

        if ( !_mesh )
            throw std::runtime_error( "Mesh is empty" );
        JeveuxVectorLong nodesRank = _mesh->getNodesRank();
        retour = nodesRank->updateValuePointer();

        if ( !_dofDescription )
            throw std::runtime_error( "Description is empty" );
        const VectorLong nodesId = _dofDescription->getNodesFromDOF();

        const int rank = getMPIRank();

        double ret = 0.0;
        for( int pos = 0; pos < taille; ++pos )
        {
            const int node_id = std::abs(nodesId[pos]);
            if((*nodesRank)[node_id-1] == rank)
                ret += ( *this )[pos] * ( *tmp )[pos];
        }

#ifdef ASTER_HAVE_MPI
        if( _mesh->isParallel() )
        {
            double ret2 = ret;
            MPIContainerUtilities utils;
            utils.all_reduce(ret2, ret, MPI_SUM);
        }
#endif

        return ret;
    }

    /**
     * @brief Size of the FieldOnNodes
     */
    ASTERINTEGER size( void ) const
    {
        return _valuesList->size();
    }


    /**
     * @brief Get DOFNumbering
     */
    BaseDOFNumberingPtr getDOFNumbering( void ) {
        return _dofNum;
    };

    /**
     * @brief Get FieldOnNodesDescription
     */
    FieldOnNodesDescriptionPtr getDescription( void ) {
        return _dofDescription;
    };


    /**
     * @brief Update field and build FieldOnNodesDescription if necessary
     */
    bool update() {
        if ( _dofNum != nullptr ) {
            _dofDescription = _dofNum->getDescription();
            _mesh = _dofNum->getMesh();
        } else if ( _dofDescription == nullptr && updateValuePointers() ) {
            typedef FieldOnNodesDescriptionClass FONDesc;
            typedef FieldOnNodesDescriptionPtr FONDescP;

            const std::string name2 = trim( ( *_reference )[1].toString() );
            _dofDescription = FONDescP( new FONDesc( name2, getMemoryType() ) );
        }
        return true;
    };

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
     */
    bool updateValuePointers() {
        bool retour = _descriptor->updateValuePointer();
        retour = ( retour && _reference->updateValuePointer() );
        retour = ( retour && _valuesList->updateValuePointer() );
        return retour;
    };

    friend class FieldBuilder;
};

template < class ValueType >
bool FieldOnNodesClass< ValueType >::printMedFile( const std::string fileName ) const {
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

/** @typedef FieldOnNodesClassReal Class d'un champ aux noeuds de double */
typedef FieldOnNodesClass< double > FieldOnNodesRealClass;

/**
 * @typedef FieldOnNodesPtrReal
 * @brief Definition d'un champ aux noeuds de double
 */
typedef boost::shared_ptr< FieldOnNodesRealClass > FieldOnNodesRealPtr;

/** @typedef FieldOnNodesLongClass Class d'une carte de long */
typedef FieldOnNodesClass< ASTERINTEGER > FieldOnNodesLongClass;

/**
 * @typedef FieldOnNodesLongPtr
 * @brief Definition d'un champ aux noeuds de long
 */
typedef boost::shared_ptr< FieldOnNodesLongClass > FieldOnNodesLongPtr;

/** @typedef FieldOnNodesClassComplex Class d'un champ aux noeuds de complexes */
typedef FieldOnNodesClass< RealComplex > FieldOnNodesComplexClass;

/**
 * @typedef FieldOnNodesComplexPtr
 * @brief Definition d'un champ aux noeuds de complexes
 */
typedef boost::shared_ptr< FieldOnNodesComplexClass > FieldOnNodesComplexPtr;

#endif /* FIELDONNODES_H_ */
