#ifndef CONTACT_PARAM_H_
#define CONTACT_PARAM_H_

/**
 * @file ContactZone.h
 * @brief Fichier entete de la class ContactZone
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

#include "Contact/ContactEnum.h"
#include "astercxx.h"

class ContactParameter {
  private:
    /** @brief Contact algorithm = ALGO_CONT */
    ContactAlgo _algo;
    /** @brief Contact algorithm = TYPE_CONT */
    ContactType _type;
    /** @brief Contact algorithm = VARIANTE */
    ContactVariant _vari;
    /** @brief Contact coefficient = COEF_CONT */
    ASTERDOUBLE _coeff;

  public:
    /**
     * @typedef ContactParameterPtr
     * @brief Pointeur intelligent vers un ContactParameter
     */
    typedef boost::shared_ptr< ContactParameter > ContactParameterPtr;

    ContactParameter()
        : _algo( ContactAlgo::Lagrangian ), _type( ContactType::Unilateral ),
          _vari( ContactVariant::Empty ), _coeff( 100. ){};

    ContactAlgo getAlgorithm() const { return _algo; };

    ContactType getType() const { return _type; };

    ContactVariant getVariant() const { return _vari; };

    ASTERDOUBLE getCoefficient() const { return _coeff; };

    void setAlgorithm( const ContactAlgo &algo ) { _algo = algo; };

    void setType( const ContactType &type ) { _type = type; };

    void setVariant( const ContactVariant &variant ) { _vari = variant; };

    void setCoefficient( const ASTERDOUBLE &coeff ) { _coeff = coeff; };
};

/**
 * @typedef ContactParameterPtr
 * @brief Pointeur intelligent vers un ContactParameter
 */
typedef boost::shared_ptr< ContactParameter > ContactParameterPtr;

class FrictionParameter {
  private:
    /** @brief Has friction ? = FROTTEMENT */
    bool _friction;
    /** @brief Friction algorithm = ALGO_FROT */
    FrictionAlgo _algo;
    /** @brief Friction algorithm = TYPE_FROT */
    FrictionType _type;
    /** @brief Friction coefficient = COEF_FROT */
    ASTERDOUBLE _coeff;

  public:
    /**
     * @typedef FrictionParameterPtr
     * @brief Pointeur intelligent vers un FrictionParameter
     */
    typedef boost::shared_ptr< FrictionParameter > FrictionParameterPtr;
};

/**
 * @typedef FrictionParameterPtr
 * @brief Pointeur intelligent vers un FrictionParameter
 */
typedef boost::shared_ptr< FrictionParameter > FrictionParameterPtr;

class PairingParameter {
  private:
  /** @brief Pairing algorithm = APPARIEMENT */
    PairingAlgo _algo;
  /** @brief Additional pairing distance = DIST_APPA */
    ASTERDOUBLE _dist_appa;

  public:
    /**
     * @typedef PairingParameterPtr
     * @brief Pointeur intelligent vers un PairingParameter
     */
    typedef boost::shared_ptr< PairingParameter > PairingParameterPtr;

    PairingParameter() : _algo(PairingAlgo::Mortar), _dist_appa(-1.0){};
};

/**
 * @typedef PairingParameterPtr
 * @brief Pointeur intelligent vers un PairingParameter
 */
typedef boost::shared_ptr< PairingParameter > PairingParameterPtr;

#endif /* CONTACT_PARAM_H_ */
