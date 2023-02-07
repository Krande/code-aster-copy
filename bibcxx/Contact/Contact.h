#ifndef CONTACT_H_
#define CONTACT_H_

/**
 * @file Contact.h
 * @brief Fichier entete de la class Contact
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

#include "astercxx.h"

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/ResultNaming.h"
#include "Utilities/CapyConvertibleValue.h"

#include <list>
#include <stdexcept>
#include <string>

class Contact : public DataStructure {
  private:
    /** @brief La SD est-elle vide ? */
    bool _isEmpty;
    /** @brief Modele */
    ModelPtr _model;
    /** @brief Ligel ".CHME.LIGRE" */
    FiniteElementDescriptorPtr _FEDesc;
    // all formulations
    /** @brief V K8 ".CHME.MODEL.NOMO" */
    JeveuxVectorChar8 _model_name;
    /** @brief V I  ".PARACI         " */
    JeveuxVectorLong _integer_params;
    /** @brief V R  ".PARACR         " */
    JeveuxVectorReal _real_params;
    /** @brief V K8 ".TYPE           " */
    JeveuxVectorChar8 _contact_type;
    // liaison_unil formulation
    /** @brief V I  ".UNILATE.NDIMCU " */
    JeveuxVectorLong _ndim_unilate;
    /** @brief V K8 ".UNILATE.CMPGCU " */
    JeveuxVectorChar8 _cmpg_unilate;
    /** @brief V K8 ".UNILATE.COED   " */
    JeveuxVectorChar8 _coed_unilate;
    /** @brief V K8 ".UNILATE.COEG   " */
    JeveuxVectorChar8 _coeg_unilate;
    /** @brief V I  ".UNILATE.LISNOE " */
    JeveuxVectorLong _lisnoe_unilate;
    /** @brief V I  ".UNILATE.POINOE " */
    JeveuxVectorLong _poinoe_unilate;
    /** @brief V R  ".UNILATE.COEFPE " */
    JeveuxVectorReal _coefpe_unilate;
    // unilateral formulation
    /** @brief  V I  ".CONTACT.NDIMCO " */
    JeveuxVectorLong _ndimco;
    // continious and discrete formulation
    /** @brief V I  ".CONTACT.METHCO " */
    JeveuxVectorLong _methco;
    /** @brief V K8 ".CONTACT.DIRAPP " */
    JeveuxVectorChar8 _dirapp;
    /** @brief V K8 ".CONTACT.DIRNOR " */
    JeveuxVectorChar8 _dirnor;
    /** @brief V K8 ".CONTACT.JFO1CO " */
    JeveuxVectorChar8 _jfo1co;
    /** @brief V K8 ".CONTACT.JFO2CO " */
    JeveuxVectorChar8 _jfo2co;
    /** @brief V R  ".CONTACT.TOLECO " */
    JeveuxVectorReal _toleco;
    /** @brief V R  ".CONTACT.JEUCOQ " */
    JeveuxVectorReal _jeucoq;
    /** @brief V R  ".CONTACT.JEUPOU " */
    JeveuxVectorReal _jeupou;

    /** @brief V I  ".CONTACT.PZONECO" */
    JeveuxVectorLong _pzoneco;
    /** @brief V I  ".CONTACT.PSUMACO" */
    JeveuxVectorLong _psumaco;
    /** @brief V I  ".CONTACT.PSUNOCO" */
    JeveuxVectorLong _psunoco;
    /** @brief V I  ".CONTACT.MAILCO " */
    JeveuxVectorLong _mailco;
    /** @brief V I  ".CONTACT.NOEUCO " */
    JeveuxVectorLong _noeuco;
    /** @brief V I  ".CONTACT.MANOCO " */
    JeveuxVectorLong _manoco;
    /** @brief V I  ".CONTACT.PMANOCO" */
    JeveuxVectorLong _pmanoco;
    /** @brief V I  ".CONTACT.NOMACO " */
    JeveuxVectorLong _nomaco;
    /** @brief V I  ".CONTACT.PNOMACO" */
    JeveuxVectorLong _pnomaco;

    /** @brief V I  ".CONTACT.PSSNOCO" */
    JeveuxVectorLong _pssnoco;
    /** @brief V I  ".CONTACT.SSNOCO " */
    JeveuxVectorLong _ssnoco;

    /** @brief V I  ".CONTACT.TYPENO " */
    JeveuxVectorLong _typeno;
    /** @brief V I  ".CONTACT.TYPEMA " */
    JeveuxVectorLong _typema;
    /** @brief V I  ".CONTACT.MAESCL " */
    JeveuxVectorLong _maescl;

    // discrete formulation
    /** @brief V R  ".CONTACT.CARADF " */
    JeveuxVectorReal _caradf;

    // continious formulation
    /** @brief V R  ".CONTACT.CARACF " */
    JeveuxVectorReal _caracf;
    /** @brief V I  ".CONTACT.PSANOFR" */
    JeveuxVectorLong _psanofr;
    /** @brief V I  ".CONTACT.SANOFR " */
    JeveuxVectorLong _sanofr;
    /** @brief V I  ".CONTACT.EXCLFR " */
    JeveuxVectorLong _exclfr;

    // xfem formulation
    /** @brief V R  ".CONTACT.CARAXF " */
    JeveuxVectorReal _caraxf;
    /** @brief V K8 ".CONTACT.XFIMAI " */
    JeveuxVectorChar8 _xfimai;
    /** @brief V K24 ".CONTACT.XNRELL " */
    JeveuxVectorChar24 _xnrell;
    /** @brief V I  ".CONTACT.MAESCX " */
    JeveuxVectorLong _maescx;

    // lac contact
    /** @brief V I  ".CONTACT.PTRDCLC" */
    JeveuxVectorLong _ptrdclc;

  public:
    /**
     * @typedef ContactPt
     * @brief Pointeur intelligent vers un Contact
     */
    typedef std::shared_ptr< Contact > ContactPtr;
    /**
     * @brief Constructeur
     */
    Contact() = delete;

    /**
     * @brief Constructeur
     */
    Contact( const std::string name, const ModelPtr model );
    /**
     * @brief Constructeur
     */
    Contact( const ModelPtr model );

    ModelPtr getModel() const;
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const;
};

/**
 * @typedef ContactPt
 * @brief Pointeur intelligent vers un Contact
 */
typedef std::shared_ptr< Contact > ContactPtr;

#endif /* CONTACT_H_ */
