#ifndef INTERSPECTRAL_H_
#define INTERSPECTRAL_H_

/**
 * @file InterspectralMatrix.h
 * @brief Fichier entete de la classe InterspectralMatrix
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "astercxx.h"

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "Supervis/ResultNaming.h"

/**
 * @class InterspectralMatrix
 * @brief Cette classe correspond a un comb_fourier
 * @author Nicolas Sellenet
 */
class InterspectralMatrix : public DataStructure {
  private:
    /** @brief Objet Jeveux '.REFE' */
    JeveuxVectorChar16 _refe;
    /** @brief Objet Jeveux '.DISC' */
    JeveuxVectorReal _disc;
    /** @brief Objet Jeveux '.VALE' */
    JeveuxCollectionReal _vale;
    /** @brief Objet Jeveux '.NUMI' */
    JeveuxVectorLong _numi;
    /** @brief Objet Jeveux '.NUMJ' */
    JeveuxVectorLong _numj;
    /** @brief Objet Jeveux '.NUME_ORDRE' */
    JeveuxVectorLong _numeOrdre;
    /** @brief Objet Jeveux '.NOEI' */
    JeveuxVectorChar8 _noei;
    /** @brief Objet Jeveux '.NOEJ' */
    JeveuxVectorChar8 _noej;
    /** @brief Objet Jeveux '.CMPI' */
    JeveuxVectorChar8 _cmpi;
    /** @brief Objet Jeveux '.CMPJ' */
    JeveuxVectorChar8 _cmpj;

    static VectorString toString( const std::vector< JeveuxChar8 > & );

  public:
    /**
     * @typedef InterspectralMatrixPtr
     * @brief Pointeur intelligent vers un InterspectralMatrix
     */
    typedef std::shared_ptr< InterspectralMatrix > InterspectralMatrixPtr;

    /**
     * @brief Constructeur
     */
    InterspectralMatrix();

    /**
     * @brief Constructeur
     */
    InterspectralMatrix( const std::string name );

    VectorLong getLineIndexes() const;
    VectorLong getColumnIndexes() const;
    VectorString getLineNodes() const;
    VectorString getColumnNodes() const;
    VectorString getLineComponents() const;
    VectorString getColumnComponents() const;
    VectorReal getNumberOfFrequencies() const;
};

/**
 * @typedef InterspectralMatrixPtr
 * @brief Pointeur intelligent vers un InterspectralMatrix
 */
typedef std::shared_ptr< InterspectralMatrix > InterspectralMatrixPtr;

#endif /* INTERSPECTRAL_H_ */
