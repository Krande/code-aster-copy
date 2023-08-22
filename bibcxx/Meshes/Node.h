/**
 * @file MeshCoordinatesField.h
 * @brief Fichier entete de la classe MeshCoordinatesField
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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

#pragma once

#include "astercxx.h"

/**
 * @class Node
 * @brief Node of a mesh
 * @author Nicolas Sellenet
 */
class Node {
  private:
    /** @brief Coordinates */
    std::array< ASTERDOUBLE, 3 > _coor;
    /** @brief Index */
    ASTERINTEGER _id;

  public:
    /**
     * @typedef MeshCoordinatesFieldPtr
     * @brief Pointeur intelligent vers un MeshCoordinatesField
     */
    typedef std::shared_ptr< Node > NodePtr;

    /**
     * @brief Constructeur
     * @param name Nom Jeveux du champ aux noeuds
     */
    Node( const ASTERINTEGER &id, const std::array< ASTERDOUBLE, 3 > &coor )
        : _id( id ), _coor( coor ) {};

    /** @brief restricted constructor (Set) and method (Get) to support pickling */
    Node( const py::tuple &tup )
        : Node( tup[0].cast< ASTERINTEGER >(), tup[1].cast< std::array< ASTERDOUBLE, 3 > >() ) {};
    py::tuple _getState() const { return py::make_tuple( _id, _coor ); };

    /**
     * @brief Get _valuesList
     */
    const std::array< ASTERDOUBLE, 3 > getValues() const { return _coor; };

    /**
     * @brief Surcharge de l'operateur []
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    ASTERDOUBLE operator[]( const ASTERINTEGER &i ) const { return _coor[i]; };

    /**
     * @brief Surcharge de l'operateur []
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    ASTERDOUBLE &operator[]( const ASTERINTEGER &i ) { return _coor[i]; };

    /**
     * @brief PlusEqual overloading
     * @return Updated field
     */
    Node &operator+=( const std::array< ASTERDOUBLE, 3 > &rhs ) {
        _coor[0] += rhs[0];
        _coor[1] += rhs[1];
        _coor[2] += rhs[2];
        return *this;
    };

    Node &operator-=( const std::array< ASTERDOUBLE, 3 > &rhs ) {
        _coor[0] -= rhs[0];
        _coor[1] -= rhs[1];
        _coor[2] -= rhs[2];
        return *this;
    };

    ASTERDOUBLE x() const { return _coor[0]; };
    ASTERDOUBLE y() const { return _coor[1]; };
    ASTERDOUBLE z() const { return _coor[2]; };

    ASTERINTEGER getId() const { return _id; };

    /**
     * @brief Size of the FieldOnNodes
     */
    ASTERINTEGER size( void ) const { return 3; }
};

/**
 * @typedef MeshCoordinatesFieldPtr
 * @brief Definition d'un champ aux noeuds de double
 */
using NodePtr = std::shared_ptr< Node >;
