#ifndef INTERSPECTRAL_H_
#define INTERSPECTRAL_H_

/**
 * @file InterspectralMatrix.h
 * @brief Fichier entete de la classe InterspectralMatrix
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

  public:
    /**
     * @typedef InterspectralMatrixPtr
     * @brief Pointeur intelligent vers un InterspectralMatrix
     */
    typedef std::shared_ptr< InterspectralMatrix > InterspectralMatrixPtr;

    /**
     * @brief Constructeur
     */
    InterspectralMatrix() : InterspectralMatrix( ResultNaming::getNewResultName() ){};

    /**
     * @brief Constructeur
     */
    InterspectralMatrix( const std::string name )
        : DataStructure( name, 8, "INTERSPECTRE" ),
          _refe( JeveuxVectorChar16( getName() + ".REFE" ) ),
          _disc( JeveuxVectorReal( getName() + ".DISC" ) ),
          _vale( JeveuxCollectionReal( getName() + ".VALE" ) ),
          _numi( JeveuxVectorLong( getName() + ".NUMI" ) ),
          _numj( JeveuxVectorLong( getName() + ".NUMJ" ) ),
          _numeOrdre( JeveuxVectorLong( getName() + ".NUME_ORDRE" ) ),
          _noei( JeveuxVectorChar8( getName() + ".NOEI" ) ),
          _noej( JeveuxVectorChar8( getName() + ".NOEJ" ) ),
          _cmpi( JeveuxVectorChar8( getName() + ".CMPI" ) ),
          _cmpj( JeveuxVectorChar8( getName() + ".CMPJ" ) ){};
};

/**
 * @typedef InterspectralMatrixPtr
 * @brief Pointeur intelligent vers un InterspectralMatrix
 */
typedef std::shared_ptr< InterspectralMatrix > InterspectralMatrixPtr;

#endif /* INTERSPECTRAL_H_ */
