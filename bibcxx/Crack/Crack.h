#ifndef CRACK_H_
#define CRACK_H_

/**
 * @file Crack.h
 * @brief Fichier entete de la classe Crack
 * @author Nicolas Pignet
 * @section LICENCE
 *   Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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

/* person_in_charge: nicolas.pignet at edf.fr */

#include "astercxx.h"

#include "DataFields/FieldOnNodes.h"
#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Supervis/ResultNaming.h"

/**
 * @class Crack
 * @brief Cette classe decrit une sd fond_fissure
 * @author Nicolas Pignet
 */
class Crack : public DataStructure {
  private:
    /** @brief Objet Jeveux '.LEVREINF.MAIL' : Vecteur (K8) contenant la liste
    des mailles de la lèvre inférieure de la fissure. */
    JeveuxVectorChar8 _levreInfMail;
    /** @brief Objet Jeveux '.NORMALE' : Vecteur de trois rééls contenant
     les composantes de la normale au plan des lèvres d'une fissure plane. */
    JeveuxVectorReal _normale;
    /** @brief Objet Jeveux '.FOND.NOEU' :  Vecteur (K8) contenant la liste
     des N noeuds ordonnés du fond de fissure. */
    JeveuxVectorChar8 _fondNoeu;
    /** @brief Objet Jeveux '.INFNORM.NOEU' : Vecteur (K8) contenant la liste
     des noeuds de la lèvre inférieure sur la direction normale au fond de fissure. */
    JeveuxVectorChar8 _infNormNoeud;
    /** @brief Objet Jeveux '.SUPNORM.NOEU' : Vecteur (K8) contenant la liste
     des noeuds de la lèvre supérieure sur la direction normale au fond de fissure. */
    JeveuxVectorChar8 _supNormNoeu;
    /** @brief Objet Jeveux '.INFNORM.NOEU2' : Vecteur (K8) contenant la liste
     des noeuds de la lèvre inférieure sur la direction normale au fond de fissure. POST_J*/
    JeveuxVectorChar8 _infNormNoeud2;
    /** @brief Objet Jeveux '.SUPNORM.NOEU2' : Vecteur (K8) contenant la liste
     des noeuds de la lèvre supérieure sur la direction normale au fond de fissure. POST_J*/
    JeveuxVectorChar8 _supNormNoeu2;
    /** @brief Objet Jeveux '.LEVRESUP.MAIL' : Vecteur (K8) contenant la liste
     des mailles de la lèvre supérieure de la fissure. */
    JeveuxVectorChar8 _levreSupMail;
    /** @brief Objet Jeveux '.INFO' : Vecteur (K24) contenant les informations sur la fissure. */
    JeveuxVectorChar24 _info;
    /** @brief Objet Jeveux '.FOND.TAILLE_R' :  Vecteur de réels contenant pour chacun
     des noeuds du fond, une estimation de la taille suivant la direction radiale,
     des mailles qui leur sont connectées. */
    JeveuxVectorReal _fondTailleR;
    /** @brief Objet Jeveux '.ABSCUR' : Vecteur de réels contenant
     les abscisses curvilignes des noeuds du fond. */
    JeveuxVectorReal _abscur;
    /** @brief Objet Jeveux '.LTNO' : Champ aux noeuds scalaire qui contient pour
     chaque noeud du maillage la valeur réelle de la level set normale à la fissure. */
    FieldOnNodesRealPtr _ltno;
    /** @brief Objet Jeveux '.LNNO' : Champ aux noeuds scalaire qui contient pour
     chaque noeud du maillage la valeur réelle de la level set tangente à la fissure. */
    FieldOnNodesRealPtr _lnno;
    /** @brief Objet Jeveux '.BASLOC' : Champ aux noeuds contenant les coordonnées
     des noeuds projetés sur le fond de fissure ainsi que les bases locales
      pour tous les noeuds du maillage.*/
    FieldOnNodesRealPtr _basLoc;
    /** @brief base local sur les noeud du fond*/
    JeveuxVectorReal _basNof;
    /** @brief Objet Jeveux '.ABSFON' : Vecteur de réels contenant
     la valeur des abscisses curvilignes des noeuds du fond de fissure. */
    JeveuxVectorReal _absfon;
    /** @brief Objet Jeveux '.COORFOND' : Vecteur de réels contenant
     les coordonnées des noeuds du fond de fissure. */
    JeveuxVectorReal _coorfond;

  public:
    /**
     * @typedef CrackPtr
     * @brief Pointeur intelligent vers un Crack
     */
    typedef std::shared_ptr< Crack > CrackPtr;

    /**
     * @brief Constructeur
     */
    Crack( const std::string name = ResultNaming::getNewResultName() );

    void updateValuePointers();
    std::string getCrackTipCellsType();
    std::string getUpperLipGroupName();
    std::string getLowerLipGroupName();
    const VectorString getCrackFrontNodes();
    const VectorString getLowerNormNodes();
    const VectorString getUpperNormNodes();
    const VectorString getLowerNormNodes2();
    const VectorString getUpperNormNodes2();
    const JeveuxVectorReal getCrackFrontBasis();
    const JeveuxVectorReal getCrackFrontPosition();
    const JeveuxVectorReal getCrackFrontAbsCurv();
    const FieldOnNodesRealPtr getCrackFrontNodeBasis();
    const JeveuxVectorReal getCrackFrontRadius();
    const JeveuxVectorReal getNormal();
    bool isSymmetric();
    std::string getConfigInit();
};

/**
 * @typedef CrackPtr
 * @brief Pointeur intelligent vers un Crack
 */
typedef std::shared_ptr< Crack > CrackPtr;

#endif /* CRACK_H_ */
