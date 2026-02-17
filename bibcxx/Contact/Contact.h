#ifndef CONTACT_H_
#define CONTACT_H_

/**
 * @file Contact.h
 * @brief Fichier entete de la class Contact
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2026  EDF www.code-aster.org
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

    /** @brief Object to model (JEVEUX name) */
    JeveuxVectorChar8 _modelName;

    /** @brief Object for type of load */
    JeveuxVectorChar8 _loadType;

    /** @brief Global parameters (integers) - PARACI" */
    JeveuxVectorLong _paraInteger;

    /** @brief Global parameters (reals) - PARACR " */
    JeveuxVectorReal _paraReal;

    // For LIAISON UNIL
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

    // General parameters for contact
    /** @brief Size and dimensions ".CONTACT.NDIMCO " */
    JeveuxVectorLong _ndimco;
    /** @brief Methods and algorithms (by zone)  ".CONTACT.METHCO " */
    JeveuxVectorLong _methco;
    /** @brief For DIRE_APPA keyword (by zone) ".CONTACT.DIRAPP " */
    JeveuxVectorChar8 _dirapp;
    /** @brief For VECT_ESCL/VECT_MAIT keywords (by zone) ".CONTACT.DIRNOR " */
    JeveuxVectorChar8 _dirnor;
    /** @brief For DIST_MAIT keyword (by zone) ".CONTACT.JFO1CO " */
    JeveuxVectorChar8 _jfo1co;
    /** @brief For DIST_ESCL keyword (by zone) ".CONTACT.JFO2CO " */
    JeveuxVectorChar8 _jfo2co;
    /** @brief For TOLE_PROJ_EXT, RESI_APPA, DIST_APPA and TOLE_INTERP keywords (by zone)
     * ".CONTACT.TOLECO " */
    JeveuxVectorReal _toleco;
    /** @brief For gap from shell/plate (by zone)  ".CONTACT.JEUCOQ " */
    JeveuxVectorReal _jeucoq;
    /** @brief For gap from beam (by zone)  ".CONTACT.JEUPOU " */
    JeveuxVectorReal _jeupou;

    /** @brief To describe mesh of surfaces (by zone) ".CONTACT.PZONECO" */
    JeveuxVectorLong _pzoneco;
    /** @brief To describe mesh of surfaces (by zone) ".CONTACT.PSUMACO" */
    JeveuxVectorLong _psumaco;
    /** @brief To describe mesh of surfaces (by zone) ".CONTACT.PSUNOCO" */
    JeveuxVectorLong _psunoco;
    /** @brief To describe mesh of surfaces (by zone) ".CONTACT.MAILCO " */
    JeveuxVectorLong _mailco;
    /** @brief To describe mesh of surfaces (by zone) ".CONTACT.NOEUCO " */
    JeveuxVectorLong _noeuco;
    /** @brief To describe mesh of surfaces (by zone) ".CONTACT.MANOCO " */
    JeveuxVectorLong _manoco;
    /** @brief To describe mesh of surfaces (by zone) ".CONTACT.PMANOCO" */
    JeveuxVectorLong _pmanoco;
    /** @brief To describe mesh of surfaces (by zone) ".CONTACT.NOMACO " */
    JeveuxVectorLong _nomaco;
    /** @brief To describe mesh of surfaces (by zone) ".CONTACT.PNOMACO" */
    JeveuxVectorLong _pnomaco;
    /** @brief To describe mesh of surfaces (by zone) ".CONTACT.PSSNOCO" */
    JeveuxVectorLong _pssnoco;
    /** @brief To describe mesh of surfaces (by zone) ".CONTACT.SSNOCO " */
    JeveuxVectorLong _ssnoco;
    /** @brief To describe mesh of surfaces (by zone) ".CONTACT.TYPENO " */
    JeveuxVectorLong _typeno;
    /** @brief To describe mesh of surfaces (by zone) ".CONTACT.TYPEMA " */
    JeveuxVectorLong _typema;
    /** @brief To describe mesh of surfaces (by zone) ".CONTACT.MAESCL " */
    JeveuxVectorLong _maescl;

    /** @brief Parameters for DISCRETE contact (by zone) ".CONTACT.CARADF " */
    JeveuxVectorReal _caradf;
    /** @brief Parameters for CONTINUE contact (by zone) ".CONTACT.CARACF " */
    JeveuxVectorReal _caracf;

    /** @brief For SANS_GROUP_NO/SANS_GROUP_NO_FR keywords (by zone) ".CONTACT.PSANOFR" */
    JeveuxVectorLong _psanofr;
    JeveuxVectorLong _sanofr;
    JeveuxVectorLong _exclfr;

    /** @brief Indirection DECOUPE_LAC/DEFI_CONTACT  ".CONTACT.PTRDCLC" */
    JeveuxVectorLong _ptrdclc;

  public:
    using ContactPtr = std::shared_ptr< Contact >;

    /** @brief No default constructor */
    Contact() = delete;

    /** @brief Constructor with given name */
    Contact( const std::string name, const ModelPtr model );

    /** @brief Constructor with automatic name */
    Contact( const ModelPtr model );

    /** @brief Get model */
    ModelPtr getModel() const;

    /** @brief Get finite element descriptor */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const;
};

using ContactPtr = std::shared_ptr< Contact >;

#endif /* CONTACT_H_ */
