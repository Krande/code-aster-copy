#ifndef CRACKTIP_H_
#define CRACKTIP_H_

/**
 * @file CrackTip.h
 * @brief Fichier entete de la classe CrackTip
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

#include "DataFields/FieldOnNodes.h"
#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Supervis/ResultNaming.h"

/**
 * @class CrackTip
 * @brief Cette classe decrit un fond_fiss
 * @author Nicolas Sellenet
 */
class CrackTip : public DataStructure {
  private:
    /** @brief Objet Jeveux '.INFO' */
    JeveuxVectorChar8 _info;
    /** @brief Objet Jeveux '.FONDFISS' */
    JeveuxVectorReal _fondFiss;
    /** @brief Objet Jeveux '.FOND.TYPE' */
    JeveuxVectorChar8 _fondType;
    /** @brief Objet Jeveux '.FOND.NOEUD' */
    JeveuxVectorChar8 _fondNoeu;
    /** @brief Objet Jeveux '.FONDINF.NOEU' */
    JeveuxVectorChar8 _fondInfNoeu;
    /** @brief Objet Jeveux '.FONDSUP.NOEU' */
    JeveuxVectorChar8 _fondSupNoeu;
    /** @brief Objet Jeveux '.FONDFISG' */
    JeveuxVectorReal _fondFisG;
    /** @brief Objet Jeveux '.NORMALE' */
    JeveuxVectorReal _normale;
    /** @brief Objet Jeveux '.BASEFOND' */
    JeveuxVectorReal _baseFond;
    /** @brief Objet Jeveux '.LTNO' */
    FieldOnNodesRealPtr _ltno;
    /** @brief Objet Jeveux '.LNNO' */
    FieldOnNodesRealPtr _lnno;
    /** @brief Objet Jeveux '.BASLOC' */
    FieldOnNodesRealPtr _basLoc;
    /** @brief Objet Jeveux '.FOND.TAILLE_R' */
    JeveuxVectorReal _fondTailleR;
    /** @brief Objet Jeveux '.DTAN_ORIGINE' */
    JeveuxVectorReal _dtanOrigine;
    /** @brief Objet Jeveux '.DTAN_EXTREMITE' */
    JeveuxVectorReal _dtanExtremite;
    /** @brief Objet Jeveux '.LEVRESUP.MAIL' */
    JeveuxVectorChar8 _levreSupMail;
    /** @brief Objet Jeveux '.SUPNORM.NOEU' */
    JeveuxVectorChar8 _supNormNoeu;
    /** @brief Objet Jeveux '.LEVREINF.MAIL' */
    JeveuxVectorChar8 _levreInfMail;
    /** @brief Objet Jeveux '.INFNORM.NOEU' */
    JeveuxVectorChar8 _infNormNoeud;

  public:
    /**
     * @typedef CrackTipPtr
     * @brief Pointeur intelligent vers un CrackTip
     */
    typedef std::shared_ptr< CrackTip > CrackTipPtr;

    /**
     * @brief Constructeur
     */
    CrackTip( const std::string name = ResultNaming::getNewResultName() );
};

/**
 * @typedef CrackTipPtr
 * @brief Pointeur intelligent vers un CrackTip
 */
typedef std::shared_ptr< CrackTip > CrackTipPtr;

#endif /* CRACKTIP_H_ */
