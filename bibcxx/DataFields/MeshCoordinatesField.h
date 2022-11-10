#ifndef MESHCOORDINATESFIELD_H_
#define MESHCOORDINATESFIELD_H_

/**
 * @file MeshCoordinatesField.h
 * @brief Fichier entete de la classe MeshCoordinatesField
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

#include "astercxx.h"

#include "DataFields/DataField.h"
#include "MemoryManager/JeveuxVector.h"

template < typename >
class FieldOnNodes;

/**
 * @class MeshCoordinatesField
 * @brief Cette classe template permet de definir un champ aux noeuds Aster
 * @author Nicolas Sellenet
 */
class MeshCoordinatesField : public DataField {
  private:
    /** @brief Vecteur Jeveux '.DESC' */
    JeveuxVectorLong _descriptor;
    /** @brief Vecteur Jeveux '.REFE' */
    JeveuxVectorChar24 _reference;
    /** @brief Vecteur Jeveux '.VALE' */
    JeveuxVectorReal _valuesList;

    using FieldOnNodesReal = FieldOnNodes< ASTERDOUBLE >;

    /**
     * @brief Return list of dof to use
     * @param local True: Use all nodes / False: Use only owned nodes
     * @param list_cmp empty: Use all cmp / keep only cmp given
     */
    VectorLong _getDOFsToUse( const VectorString &list_cmp ) const;

  public:
    /**
     * @typedef MeshCoordinatesFieldPtr
     * @brief Pointeur intelligent vers un MeshCoordinatesField
     */
    typedef std::shared_ptr< MeshCoordinatesField > MeshCoordinatesFieldPtr;

    /**
     * @brief Constructeur
     * @param name Nom Jeveux du champ aux noeuds
     */
    MeshCoordinatesField( const std::string &name )
        : DataField( name, "CHAM_NO" ),
          _descriptor( JeveuxVectorLong( getName() + ".DESC" ) ),
          _reference( JeveuxVectorChar24( getName() + ".REFE" ) ),
          _valuesList( JeveuxVectorReal( getName() + ".VALE" ) ){};

    /**
     * @brief copy constructeur
     * @param coordField MeshCoordinatesField
     */
    MeshCoordinatesField( const MeshCoordinatesField &coordField )
        : MeshCoordinatesField( ResultNaming::getNewResultName() ) {
        *( _descriptor ) = *( coordField._descriptor );
        *( _reference ) = *( coordField._reference );
        *( _valuesList ) = *( coordField._valuesList );
    };

    bool exists() const { return _valuesList->exists(); };

    /**
     * @brief Shorthand + operator assignement
     * @return Updated field
     */
    MeshCoordinatesField &operator+=( const FieldOnNodesReal &rhs );

    MeshCoordinatesField &operator+=( const MeshCoordinatesField &rhs ) {
        ( *_valuesList ) += ( *rhs._valuesList );
        return *this;
    };

    friend MeshCoordinatesField operator+( MeshCoordinatesField lhs, const FieldOnNodesReal &rhs ) {
        lhs += rhs;
        return lhs;
    };

    friend MeshCoordinatesField operator+( const FieldOnNodesReal &lhs,
                                           const MeshCoordinatesField &rhs ) {
        return rhs + lhs;
    };

    friend MeshCoordinatesField operator+( MeshCoordinatesField lhs,
                                           const MeshCoordinatesField &rhs ) {
        lhs += rhs;
        return lhs;
    };

    MeshCoordinatesField &operator-=( const FieldOnNodesReal &rhs );

    MeshCoordinatesField &operator-=( const MeshCoordinatesField &rhs ) {
        ( *_valuesList ) -= ( *rhs._valuesList );
        return *this;
    };

    friend MeshCoordinatesField operator-( MeshCoordinatesField lhs, const FieldOnNodesReal &rhs ) {
        lhs -= rhs;
        return lhs;
    };

    friend MeshCoordinatesField operator-( const FieldOnNodesReal &lhs,
                                           const MeshCoordinatesField &rhs ) {
        return -( rhs - lhs );
    };

    friend MeshCoordinatesField operator-( MeshCoordinatesField lhs,
                                           const MeshCoordinatesField &rhs ) {
        lhs -= rhs;
        return lhs;
    };

    /**
     * @brief TimesEqual overloading
     * @return Updated field
     */
    MeshCoordinatesField &operator*=( const ASTERDOUBLE &scal ) {

        ( *_valuesList ) *= scal;

        return *this;
    };

    /**
     * @brief Unary Minus overloading
     * @return Updated field
     */
    MeshCoordinatesField operator-() const {
        MeshCoordinatesField tmp( *this );
        ( *tmp._valuesList ) *= double( -1.0 );
        return tmp;
    };

    /**
     * @brief Multiply by a scalar on right overloading
     * @return New field
     */
    friend MeshCoordinatesField operator*( MeshCoordinatesField lhs, const ASTERDOUBLE &scal ) {

        lhs *= scal;
        return lhs;
    };

    friend MeshCoordinatesField operator*( const ASTERDOUBLE &scal,
                                           const MeshCoordinatesField &rhs ) {
        return rhs * scal;
    };

    /**
     * @brief Assignement operator
     * @param coordField MeshCoordinatesField
     */
    MeshCoordinatesField &operator=( const MeshCoordinatesField &coordField ) {
        *( _descriptor ) = *( coordField._descriptor );
        *( _reference ) = *( coordField._reference );
        *( _valuesList ) = *( coordField._valuesList );
        return *this;
    };

    /**
     * @brief Get _descriptor
     */
    const JeveuxVectorLong getDescriptor() const { return _descriptor; };

    /**
     * @brief Get _reference
     */
    const JeveuxVectorChar24 getReference() const { return _reference; };

    /**
     * @brief Get _valuesList
     */
    const JeveuxVectorReal getValues() const { return _valuesList; };

    /**
     * @brief Surcharge de l'operateur []
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    ASTERDOUBLE operator[]( const ASTERINTEGER &i ) const { return _valuesList->operator[]( i ); };

    /**
     * @brief Size of the FieldOnNodes
     */
    ASTERINTEGER size( void ) const { return _valuesList->size(); }

    /**
     * @brief duplicate
     * @return MeshCoordinatesField
     */
    MeshCoordinatesField duplicate() { return *this; };

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
     */
    void updateValuePointers() const {
        _descriptor->updateValuePointer();
        _reference->updateValuePointer();
        _valuesList->updateValuePointer();
    };
};

/**
 * @typedef MeshCoordinatesFieldPtr
 * @brief Definition d'un champ aux noeuds de double
 */
using MeshCoordinatesFieldPtr = std::shared_ptr< MeshCoordinatesField >;

/**
 * @typedef MeshCoordinatesFieldPtr
 * @brief Definition d'un champ aux noeuds de double
 */
using ConstMeshCoordinatesFieldPtr = std::shared_ptr< const MeshCoordinatesField >;

#endif /* MESHCOORDINATESFIELD_H_ */
