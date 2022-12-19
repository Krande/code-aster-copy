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

#include "aster_fort_petsc.h"
#include "aster_fort_superv.h"
#include "aster_fort_utils.h"
#include "astercxx.h"

#include "DataFields/DataField.h"
#include "DataFields/MeshCoordinatesField.h"
#include "DataFields/SimpleFieldOnNodes.h"
#include "MemoryManager/JeveuxAllowedTypes.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/BaseMesh.h"
#include "Numbering/DOFNumbering.h"
#include "Numbering/FieldOnNodesDescription.h"
#include "Numbering/ParallelDOFNumbering.h"
#include "ParallelUtilities/AsterMPI.h"
#include "PythonBindings/LogicalUnitManager.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/Exceptions.h"
#include "Utilities/Blas.h"
#ifdef ASTER_HAVE_PETSC
#include <petscvec.h>
#endif

#include <typeinfo>

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
    typedef std::shared_ptr< SimpleFieldOnNodesValueType > SimpleFieldOnNodesValueTypePtr;

    /** @brief Vecteur Jeveux '.DESC' */
    JeveuxVectorLong _descriptor;
    /** @brief Vecteur Jeveux '.REFE' */
    JeveuxVectorChar24 _reference;
    /** @brief Vecteur Jeveux '.VALE' */
    JeveuxVector< ValueType > _values;
    /** @brief Dof description */
    FieldOnNodesDescriptionPtr _dofDescription;
    /** @brief Support mesh */
    BaseMeshPtr _mesh;

    /**
     * @brief Return list of dof to use
     * @param sameRank True: Use only owned nodes / False: Use all nodes
     * @param list_cmp empty: Use all cmp / keep only cmp given
     */
    VectorLong _getDOFsToUse( const bool sameRank, const VectorString &list_cmp ) const {
        const bool all_nodes = !sameRank;
        const auto rank = getMPIRank();
        if ( !_mesh )
            raiseAsterError( "Mesh is empty" );
        const JeveuxVectorLong nodesRank = _mesh->getNodesRank();
        nodesRank->updateValuePointer();

        const bool all_cmp = list_cmp.empty();
        std::map< ASTERINTEGER, std::string > map_cmp;
        if ( !all_cmp ) {
            auto name2num = this->getComponentsName2Number();
            for ( auto &name : list_cmp ) {
                auto name_trim = trim( name );
                if ( name2num.count( name_trim ) > 0 ) {
                    map_cmp[name2num[name_trim]] = name_trim;
                }
            }
        }

        const auto descr = this->getNodesAndComponentsNumberFromDOF();

        auto taille = this->size();
        VectorLong dofUsed;
        dofUsed.reserve( taille );

        for ( auto dof = 0; dof < taille; ++dof ) {
            const auto node_id = std::abs( descr[dof].first );
            bool l_keep_cmp = all_cmp || map_cmp.count( descr[dof].second ) > 0;
            bool l_keep_node = all_nodes || ( *nodesRank )[node_id] == rank;
            if ( l_keep_node && l_keep_cmp )
                dofUsed.push_back( dof );
        }

        return dofUsed;
    }

  public:
    /** @typedef FieldOnNodesPtr */
    typedef std::shared_ptr< FieldOnNodes > FieldOnNodesPtr;

    /**
     * @brief Constructor
     * @param name Jeveux name of the field
     */
    FieldOnNodes( const std::string name )
        : DataField( name, "CHAM_NO" ),
          _descriptor( JeveuxVectorLong( getName() + ".DESC" ) ),
          _reference( JeveuxVectorChar24( getName() + ".REFE" ) ),
          _values( JeveuxVector< ValueType >( getName() + ".VALE" ) ),
          _dofDescription( nullptr ),
          _mesh( nullptr ){};

    /** @brief Constructor with automatic name */
    FieldOnNodes() : FieldOnNodes( DataStructureNaming::getNewName() ){};

    /** @brief Copy constructor */
    FieldOnNodes( const std::string &name, const FieldOnNodes &toCopy ) : FieldOnNodes( name ) {
        // JeveuxVector to be duplicated
        *( _descriptor ) = *( toCopy._descriptor );
        *( _reference ) = *( toCopy._reference );
        *( _values ) = *( toCopy._values );
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
        _values = other._values;
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
    FieldOnNodes( const FiniteElementDescriptorPtr &fed, const std::string &localMode )
        : FieldOnNodes() {

        if ( !fed ) {
            raiseAsterError( " FEDescriptor is null" );
        }

        // Create numbering from Local mode (see Cata)
        auto dofNume = std::make_shared< DOFNumbering >();

        dofNume->addFiniteElementDescriptor( fed );
        dofNume->computeNumberingWithLocalMode( localMode );

        _dofDescription = dofNume->getDescription();
        _mesh = dofNume->getMesh();

        const auto intType = AllowedFieldType< ValueType >::numTypeJeveux;
        CALLO_VTCREB_WRAP( getName(), JeveuxMemoryTypesNames[Permanent], JeveuxTypesNames[intType],
                           dofNume->getName() );
    };

    /**
     * @brief Constructor with DOFNumbering
     */
    FieldOnNodes( const ModelPtr &model ) : FieldOnNodes() {

        if ( !model ) {
            raiseAsterError( " Model is null" );
        }

        // Create numbering from Local mode (see Cata)
        auto dofNume = std::make_shared< DOFNumbering >();

        dofNume->setModel( model );
        dofNume->computeNumbering();

        _dofDescription = dofNume->getDescription();
        _mesh = dofNume->getMesh();

        const auto intType = AllowedFieldType< ValueType >::numTypeJeveux;
        CALLO_VTCREB_WRAP( getName(), JeveuxMemoryTypesNames[Permanent], JeveuxTypesNames[intType],
                           dofNume->getName() );
    };

    /**
     * @brief Constructor with DOFNumbering
     */
    FieldOnNodes( const BaseDOFNumberingPtr &dofNum ) : FieldOnNodes() {

        if ( !dofNum ) {
            raiseAsterError( " DOFNumbering is null" );
        }

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
          _values( toCopy->getValues() ),
          _dofDescription( nullptr ),
          _mesh( nullptr ){};

    /**
     * @brief Surcharge de l'operateur []
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    ValueType &operator[]( ASTERINTEGER i ) { return _values->operator[]( i ); };

    /**
     * @brief Surcharge de l'operateur []
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    const ValueType &operator[]( ASTERINTEGER i ) const {
        return const_cast< ValueType & >(
            const_cast< FieldOnNodes< ValueType > * >( this )->operator[]( i ) );
    };

    /**
     * @brief Check if fields are OK for +, +=, ...
     * @return true if compatible
     */
    bool isSimilarTo( const FieldOnNodes< ValueType > &tmp2 ) {
        CALL_JEMARQ();
        bool similar = ( ( *_descriptor ) == ( *tmp2._descriptor ) );
        similar = ( similar && ( this->_reference->size() == tmp2._reference->size() ) );
        similar = ( similar && ( this->_values->size() == tmp2._values->size() ) );
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
     * @brief Scale a vector by a diagonal matrix stored as a vector
     * @return Updated field
     */
    FieldOnNodes< ValueType > &scale( const VectorReal &vect ) {
        _values->updateValuePointer();
        if ( _values->size() != vect.size() )
            raiseAsterError( "Field and vector have incompatible shapes" );
        for ( std::size_t i = 0; i < vect.size(); ++i ) {
            ( *_values )[i] = ( *_values )[i] * vect[i];
        }
        return *this;
    };

    /**
     * @brief PlusEqual overloading
     * @return Updated field
     */
    FieldOnNodes< ValueType > &operator+=( FieldOnNodes< ValueType > const &rhs ) {
        if ( !this->isSimilarTo( rhs ) )
            raiseAsterError( "Fields have incompatible shapes" );
        ( *_values ) += ( *rhs._values );
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
        ( *_values ) -= ( *rhs._values );
        return *this;
    };

    /**
     * @brief TimesEqual overloading
     * @return Updated field
     */
    FieldOnNodes< ValueType > &operator*=( const ASTERDOUBLE &scal ) {
        // AsterBLAS::scal(scal, _values);

        ( *_values ) *= scal;

        return *this;
    };

    /**
     * @brief Unary Minus overloading
     * @return Updated field
     */
    FieldOnNodes< ValueType > operator-() const {
        FieldOnNodes< ValueType > tmp( *this );
        ( *tmp._values ) *= ValueType( -1 );
        return tmp;
    };

    // NOTE ON FOLLOWING OPERATOR OVERLOADING
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

        lhs *= scal;
        return lhs;
    };

    /**
     * @brief Multiply by a scalar on left overloading
     * @return New field
     */
    friend FieldOnNodes< ValueType > operator*( const ASTERDOUBLE &scal,
                                                const FieldOnNodes< ValueType > &rhs ) {
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
     * @brief Returns a vector with node index and component name for each DOFs
     */
    std::vector< std::pair< ASTERINTEGER, std::string > >
    getNodesAndComponentsFromDOF( const bool local = true ) const {
        if ( !_dofDescription )
            raiseAsterError( "Description is empty" );

        auto nodeAndComponentsNumberFromDOF = _dofDescription->getNodeAndComponentsNumberFromDOF();
        nodeAndComponentsNumberFromDOF->updateValuePointer();

        const ASTERINTEGER nb_eq = _dofDescription->getNumberOfDofs();

        auto num2name = this->getComponentsNumber2Name();

        std::vector< std::pair< ASTERINTEGER, std::string > > ret;
        ret.reserve( nb_eq );

        for ( ASTERINTEGER i_eq = 0; i_eq < nb_eq; i_eq++ ) {
            auto node_id = ( *nodeAndComponentsNumberFromDOF )[2 * i_eq] - 1;
            auto cmp = abs( ( *nodeAndComponentsNumberFromDOF )[2 * i_eq + 1] );
            std::string cmp_name = " ";
            if ( cmp > 0 ) {
                cmp_name = num2name[cmp];
            }
            ret.push_back( std::make_pair( node_id, cmp_name ) );
        }

#ifdef ASTER_HAVE_MPI
        if ( !local && _mesh->isParallel() ) {

            auto mapLG = _mesh->getLocalToGlobalMapping();
            mapLG->updateValuePointer();
            for ( ASTERINTEGER i_eq = 0; i_eq < nb_eq; i_eq++ ) {
                auto node_id = ret[i_eq].first;
                ret[i_eq].first = ( *mapLG )[node_id];
            }
        }
#endif

        return ret;
    };

    /**
     * @brief Returns a vector with node index and component name for each DOFs
     */
    std::vector< std::pair< ASTERINTEGER, ASTERINTEGER > >
    getNodesAndComponentsNumberFromDOF( const bool local = true ) const {
        if ( !_dofDescription )
            raiseAsterError( "Description is empty" );

        auto nodeAndComponentsNumberFromDOF = _dofDescription->getNodeAndComponentsNumberFromDOF();
        nodeAndComponentsNumberFromDOF->updateValuePointer();

        const ASTERINTEGER nb_eq = _dofDescription->getNumberOfDofs();

        std::vector< std::pair< ASTERINTEGER, ASTERINTEGER > > ret;
        ret.reserve( nb_eq );

        for ( int i_eq = 0; i_eq < nb_eq; i_eq++ ) {
            auto node_id = ( *nodeAndComponentsNumberFromDOF )[2 * i_eq] - 1;
            auto cmp = ( *nodeAndComponentsNumberFromDOF )[2 * i_eq + 1];
            ret.push_back( std::make_pair( node_id, cmp ) );
        }

#ifdef ASTER_HAVE_MPI
        if ( !local && _mesh->isParallel() ) {

            auto mapLG = _mesh->getLocalToGlobalMapping();
            mapLG->updateValuePointer();
            for ( ASTERINTEGER i_eq = 0; i_eq < nb_eq; i_eq++ ) {
                auto node_id = ret[i_eq].first;
                ret[i_eq].first = ( *mapLG )[node_id];
            }
        }
#endif

        return ret;
    };

    /**
     * @brief Maps between name of components and the nimber
     */
    std::map< std::string, ASTERINTEGER > getComponentsName2Number() const {
        std::map< std::string, ASTERINTEGER > ret;

        JeveuxVectorChar8 cmp_name( "&CMP_NAME" );
        JeveuxVectorLong cmp( "&CMP" );
        ASTERINTEGER ncmp;

        CALL_UTNCMP2( getName().c_str(), &ncmp, cmp->getName().c_str(),
                      cmp_name->getName().c_str() );

        cmp_name->updateValuePointer();
        cmp->updateValuePointer();

        for ( int icmp = 0; icmp < ncmp; icmp++ ) {
            ret[trim( ( *cmp_name )[icmp] )] = ( *cmp )[icmp];
        }

        return ret;
    };

    std::map< ASTERINTEGER, std::string > getComponentsNumber2Name() const {
        std::map< ASTERINTEGER, std::string > ret;

        auto name2number = this->getComponentsName2Number();
        for ( auto &[name, num] : name2number ) {
            ret[num] = name;
        }

        return ret;
    };

    VectorString getComponents() const {

        JeveuxVectorChar8 cmp( "&CMP" );
        ASTERINTEGER ncmp;

        CALL_UTNCMP( getName().c_str(), &ncmp, cmp->getName().c_str() );

        cmp->updateValuePointer();

        VectorString cmps;
        cmps.reserve( cmp->size() );

        for ( auto &cm : cmp ) {
            cmps.push_back( trim( cm.toString() ) );
        }

        return cmps;
    }

    ASTERINTEGER getNumberOfComponents() const { return getComponents().size(); }

    /**
     * @brief Allouer un champ au noeud à partir d'un autre
     * @return renvoit true
     */
    bool allocateFrom( const FieldOnNodes< ValueType > &tmp ) {
        this->_descriptor->deallocate();
        this->_reference->deallocate();
        this->_values->deallocate();

        this->_descriptor->allocate( tmp._descriptor->size() );
        this->_reference->allocate( tmp._reference->size() );
        this->_values->allocate( tmp._values->size() );
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
        toReturn->build();
        return toReturn;
    };

    bool exists() const {
        return _reference->exists() && _descriptor->exists() && _values->exists();
    };

    /**
     * @brief Get mesh
     */
    BaseMeshPtr getMesh() const { return _mesh; };

    bool printMedFile( const std::string fileName, bool local = true ) const;

    /**
     * @brief Import a PETSc vector into a ParallelFieldOnNodes
     *
     * @param dofNmbrg The numbering of the DOFs
     * @param vec The PETSc vector
     * @param scaling The scaling of the Lagrange DOFs
     */
#ifdef ASTER_HAVE_PETSC
    void fromPetsc( const DOFNumbering &dofNmbrg, const Vec &vec, const ASTERDOUBLE &scaling ) {
        if ( get_sh_jeveux_status() == 1 ) {
            CALLO_VECT_ASSE_FROM_PETSC( getName(), dofNmbrg.getName(), &vec, &scaling );
            _values->updateValuePointer();
        }
    };
    void fromPetsc( const ParallelDOFNumbering &dofNmbrg, const Vec &vec,
                    const ASTERDOUBLE &scaling ) {
        if ( get_sh_jeveux_status() == 1 ) {
            CALLO_VECT_ASSE_FROM_PETSC( getName(), dofNmbrg.getName(), &vec, &scaling );
            _values->updateValuePointer();
        }
    };
    void fromPetsc( const DOFNumbering &dofNmbrg, const Vec &vec ) {
        if ( get_sh_jeveux_status() == 1 ) {
            const auto dummy_scaling = 1.;
            CALLO_VECT_ASSE_FROM_PETSC( getName(), dofNmbrg.getName(), &vec, &dummy_scaling );
            _values->updateValuePointer();
        }
    };
    void fromPetsc( const ParallelDOFNumbering &dofNmbrg, const Vec &vec ) {
        if ( get_sh_jeveux_status() == 1 ) {
            const auto dummy_scaling = 1.;
            CALLO_VECT_ASSE_FROM_PETSC( getName(), dofNmbrg.getName(), &vec, &dummy_scaling );
            _values->updateValuePointer();
        }
    };
    void fromPetsc( const Vec &vec, const ASTERDOUBLE &scaling ) {
        if ( get_sh_jeveux_status() == 1 ) {
            if ( getMesh()->isParallel() ) {
                raiseAsterError( "dofNumbering must be provided" );
            }
            const std::string dummy_nbg = " ";
            CALLO_VECT_ASSE_FROM_PETSC( getName(), dummy_nbg, &vec, &scaling );
            _values->updateValuePointer();
        }
    };
    void fromPetsc( const Vec &vec ) {
        if ( get_sh_jeveux_status() == 1 ) {
            if ( getMesh()->isParallel() ) {
                raiseAsterError( "dofNumbering must be provided" );
            }
            const auto dummy_scaling = 1.;
            const std::string dummy_nbg = " ";
            CALLO_VECT_ASSE_FROM_PETSC( getName(), dummy_nbg, &vec, &dummy_scaling );
            _values->updateValuePointer();
        }
    };
#endif

    /**
     * @brief Set the Values object
     *
     * @param value Value to affect
     */
    void setValues( const ValueType &value ) {
        CALL_JEMARQ();
        _values->updateValuePointer();
        _values->assign( value );
        CALL_JEDEMA();
    };

    void setValues( const std::vector< ValueType > &values ) {
        AS_ASSERT( values.size() == size() );

        *_values = values;
    };

    void applyLagrangeScaling( const ValueType scaling ) {
        _values->updateValuePointer();

        const auto descr = this->getNodesAndComponentsNumberFromDOF();

        auto nbDofs = this->size();

        for ( ASTERINTEGER dof = 0; dof < nbDofs; dof++ ) {
            if ( descr[dof].second < 0 ) {
                ( *this )[dof] = ( *this )[dof] * scaling;
            }
        }
    };

    void setValues( const std::map< std::string, ValueType > &values ) {
        _values->updateValuePointer();

        VectorString list_cmp;
        for ( auto &[key, val] : values ) {
            list_cmp.push_back( trim( key ) );
        }

        auto num2name = this->getComponentsNumber2Name();
        const auto descr = this->getNodesAndComponentsNumberFromDOF();

        auto nbDofs = this->size();

        for ( ASTERINTEGER dof = 0; dof < nbDofs; dof++ ) {
            auto search = values.find( num2name[descr[dof].second] );
            if ( search != values.end() ) {
                ( *this )[dof] = search->second;
            }
        }
    };

    /**
     * @brief Get values of the field
     *
     */
    const JeveuxVector< ValueType > &getValues() const { return _values; }

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
     * @brief Compute norm
     * @param normType Type of norm ("NORM_1","NORM_2","NORM_INFINITY")
     */
    ASTERDOUBLE norm( const std::string normType, VectorString list_cmp = VectorString() ) const {

        if constexpr ( !std::is_same_v< ValueType, ASTERDOUBLE > &&
                       !std::is_same_v< ValueType, ASTERINTEGER > &&
                       !std::is_same_v< ValueType, ASTERCOMPLEX > ) {
            raiseAsterError( " norm method not defined for type " +
                             std::string( typeid( ValueType ).name() ) );
        }
        CALL_JEMARQ();
        ASTERDOUBLE norme = 0.0;
        _values->updateValuePointer();
        auto dofUsed = this->_getDOFsToUse( true, list_cmp );

        if ( normType == "NORM_1" ) {
            for ( auto &dof : dofUsed ) {
                norme += std::abs( ( *this )[dof] );
            }
        } else if ( normType == "NORM_2" ) {
            for ( auto &dof : dofUsed ) {
                norme += std::pow( std::abs( ( *this )[dof] ), 2 );
            }
        } else if ( normType == "NORM_INFINITY" ) {
            for ( auto &dof : dofUsed ) {
                norme = std::max( norme, std::abs( ( *this )[dof] ) );
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
    ValueType dot( const FieldOnNodesPtr &tmp ) const {
        CALL_JEMARQ();
        tmp->updateValuePointers();
        _values->updateValuePointer();
        const auto taille = _values->size();

        if ( taille != tmp->size() )
            raiseAsterError( "Incompatible size" );

        auto dofUsed = this->_getDOFsToUse( true, VectorString() );

        ValueType ret;
        if constexpr ( std::is_same_v< ValueType, ASTERDOUBLE > ||
                       std::is_same_v< ValueType, ASTERINTEGER > ) {
            ret = (ValueType)0.0;
            for ( auto &dof : dofUsed ) {
                ret += ( *this )[dof] * ( *tmp )[dof];
            }
        } else if constexpr ( std::is_same_v< ValueType, ASTERCOMPLEX > ) {
            ret = (ValueType)0.0;
            for ( auto &dof : dofUsed ) {
                ret += ( *this )[dof] * std::conj( ( *tmp )[dof] );
            }
        } else {
            raiseAsterError( " dot method not defined for type " +
                             std::string( typeid( ValueType ).name() ) );
        }

#ifdef ASTER_HAVE_MPI
        if ( _mesh->isParallel() ) {
            ValueType ret2 = ret;
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
     * @brief Get physical quantity
     */
    std::string getPhysicalQuantity() const {
        const std::string typeco( "CHAM_NO" );
        ASTERINTEGER repi = 0, ier = 0;
        JeveuxChar32 repk( " " );
        const std::string arret( "F" );
        const std::string questi( "NOM_GD" );
        CALLO_DISMOI( questi, this->getName(), typeco, &repi, repk, arret, &ier );
        auto retour = trim( repk.toString() );
        return retour;
    }

    /**
     * @brief Size of the FieldOnNodes
     */
    ASTERINTEGER size( void ) const { return _values->size(); }

    /**
     * @brief Get FieldOnNodesDescription
     */
    FieldOnNodesDescriptionPtr getDescription( void ) const { return _dofDescription; };

    /**
     * @brief Update field and build FieldOnNodesDescription if necessary
     */
    bool build() {
        if ( !_dofDescription ) {
            CALL_JEMARQ();

            _reference->updateValuePointer();
            const std::string name2 = trim( ( *_reference )[1].toString() );
            if ( !name2.empty() ) {
                _dofDescription = std::make_shared< FieldOnNodesDescription >( name2 );
            }
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
        _values->updateValuePointer();
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

    ASTERINTEGER op = 39;
    CALL_EXECOP( &op );

    return true;
};

/** @typedef FieldOnNodesReal */
using FieldOnNodesReal = FieldOnNodes< ASTERDOUBLE >;
using FieldOnNodesRealPtr = std::shared_ptr< FieldOnNodesReal >;

/** @typedef FieldOnNodesLong */
using FieldOnNodesLong = FieldOnNodes< ASTERINTEGER >;
using FieldOnNodesLongPtr = std::shared_ptr< FieldOnNodesLong >;

/** @typedef FieldOnNodesComplex*/
using FieldOnNodesComplex = FieldOnNodes< ASTERCOMPLEX >;
using FieldOnNodesComplexPtr = std::shared_ptr< FieldOnNodesComplex >;

/** @typedef FieldOnNodesChar8 */
using FieldOnNodesChar8 = FieldOnNodes< JeveuxChar8 >;
using FieldOnNodesChar8Ptr = std::shared_ptr< FieldOnNodesChar8 >;

#endif /* FIELDONNODES_H_ */
