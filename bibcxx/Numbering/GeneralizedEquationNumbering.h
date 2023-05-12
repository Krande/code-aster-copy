/**
 * @file GeneralizedEquationNumbering.h
 * @brief Fichier entete de la classe GeneralizedEquationNumbering
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

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"

/**
 * @class GeneralizedEquationNumbering
 * @brief This class describes the structure of dof stored in a field on nodes
 * @author Nicolas Sellenet
 */
class GeneralizedEquationNumbering : public DataStructure {
    /** @brief Objet Jeveux '.DESC' */
    JeveuxVectorLong _desc;
    /** @brief Objet Jeveux '.NEQU' */
    JeveuxVectorLong _nequ;
    /** @brief Objet Jeveux '.REFN' */
    JeveuxVectorChar24 _refn;
    /** @brief Objet Jeveux '.DELG' */
    JeveuxVectorLong _delg;
    /** @brief Objet Jeveux '.ORIG' */
    JeveuxCollectionLong _orig;
    /** @brief Objet Jeveux '.PRNO' */
    JeveuxCollectionLong _componentsOnNodes;
    /** @brief Objet Jeveux '.LILI' */
    NamesMapChar24 _namesOfGroupOfCells;
    /** @brief Objet Jeveux '.NUEQ' */
    JeveuxVectorLong _indexationVector;
    /** @brief Objet Jeveux '.DEEQ' */
    JeveuxVectorLong _nodeAndComponentsIdFromDOF;

  public:
    /**
     * @brief Constructeur
     * @param name nom souhaité de la sd (utile pour le GeneralizedEquationNumbering
     * d'une sd_resu)
     */
    GeneralizedEquationNumbering( const std::string name );

    /**
     * @brief Constructeur
     */
    GeneralizedEquationNumbering();

    /**
     * @brief Destructor
     */
    ~GeneralizedEquationNumbering() {};
};

typedef std::shared_ptr< GeneralizedEquationNumbering > GeneralizedEquationNumberingPtr;
