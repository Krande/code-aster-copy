#ifndef CONTACT_PARAM_H_
#define CONTACT_PARAM_H_

/**
 * @file ContactZone.h
 * @brief Fichier entete de la class ContactZone
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

#include "astercxx.h"

#include "aster_pybind.h"

#include "Contact/ContactEnum.h"
#include "Discretization/ElementaryCharacteristics.h"
#include "Loads/ListOfLoads.h"

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
    /** @brief Jacobian computation = TYPE_MATR_TANG */
    JacobianType _jacType;

  public:
    /**
     * @typedef ContactParameterPtr
     * @brief Pointeur intelligent vers un ContactParameter
     */
    typedef std::shared_ptr< ContactParameter > ContactParameterPtr;

    /**
     * @brief Constructeur
     */
    ContactParameter()
        : _algo( ContactAlgo::Lagrangian ),
          _type( ContactType::Unilateral ),
          _vari( ContactVariant::Empty ),
          _coeff( 100. ),
          _jacType( JacobianType::Analytical ) {};

    /** @brief restricted constructor (Set) and method (Get) to support pickling */
    ContactParameter( const py::tuple &tup );
    py::tuple _getState() const;

    ContactAlgo getAlgorithm() const { return _algo; };

    ContactType getType() const { return _type; };

    ContactVariant getVariant() const { return _vari; };

    ASTERDOUBLE getCoefficient() const { return _coeff; };

    JacobianType getJacobianType() const { return _jacType; };

    void setAlgorithm( const ContactAlgo &algo ) { _algo = algo; };

    void setType( const ContactType &type ) { _type = type; };

    void setVariant( const ContactVariant &variant ) { _vari = variant; };

    void setCoefficient( const ASTERDOUBLE &coeff ) { _coeff = coeff; };

    void setJacobianType( const JacobianType &type ) { _jacType = type; };
};

/**
 * @typedef ContactParameterPtr
 * @brief Pointeur intelligent vers un ContactParameter
 */
typedef std::shared_ptr< ContactParameter > ContactParameterPtr;

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
    /** @brief TRESCA coefficient = TRESCA */
    ASTERDOUBLE _tresca;
    /** @brief COULOMB coefficient = COULOMB */
    ASTERDOUBLE _coulomb;

  public:
    /**
     * @typedef FrictionParameterPtr
     * @brief Pointeur intelligent vers un FrictionParameter
     */
    typedef std::shared_ptr< FrictionParameter > FrictionParameterPtr;

    /**
     * @brief Constructeur
     */

    FrictionParameter()
        : _friction( false ),
          _algo( FrictionAlgo::Lagrangian ),
          _type( FrictionType::Without ),
          _coeff( 100. ),
          _tresca( -1. ),
          _coulomb( -1. ) {};

    /** @brief restricted constructor (Set) and method (Get) to support pickling */
    FrictionParameter( const py::tuple &tup );
    py::tuple _getState() const;

    FrictionAlgo getAlgorithm() const { return _algo; };

    FrictionType getType() const { return _type; };

    ASTERDOUBLE getCoefficient() const { return _coeff; };

    ASTERDOUBLE getTresca() const { return _tresca; };

    ASTERDOUBLE getCoulomb() const { return _coulomb; };

    void setAlgorithm( const FrictionAlgo &algo ) { _algo = algo; };

    void setType( const FrictionType &type ) { _type = type; };

    void setCoefficient( const ASTERDOUBLE &coeff ) { _coeff = coeff; };

    void setTresca( const ASTERDOUBLE &tresca ) { _tresca = tresca; };

    void setCoulomb( const ASTERDOUBLE &coulomb ) { _coulomb = coulomb; };

    void enableFriction( const bool &friction ) { _friction = friction; }

    bool hasFriction() const { return _friction; }
};

/**
 * @typedef FrictionParameterPtr
 * @brief Pointeur intelligent vers un FrictionParameter
 */
typedef std::shared_ptr< FrictionParameter > FrictionParameterPtr;

class PairingParameter {
  private:
    /** @brief Pairing algorithm = APPARIEMENT */
    PairingAlgo _algo;
    /** @brief Additional pairing distance = COEF_MULT_APPA */
    ASTERDOUBLE _dist_ratio;
    /** @brief initial contact state = CONTACT_INIT */
    InitialState _cont_init;

    /** @brief fictive distance function = DIST_SUPP */
    GenericFunctionPtr _dist_supp;
    /** @brief if fictive distance for beam = DIST_POUTRE */
    bool _beam;
    /** @brief if fictive distance for shell = DIST_COQUE */
    bool _shell;
    /** @brief structural element characteristics = CARA_ELEM */
    ElementaryCharacteristicsPtr _cara;
    /** @brief structural element characteristics = DIST_SUPP */

  public:
    /**
     * @typedef PairingParameterPtr
     * @brief Pointeur intelligent vers un PairingParameter
     */
    typedef std::shared_ptr< PairingParameter > PairingParameterPtr;

    /**
     * @brief Constructeur
     */
    PairingParameter()
        : _algo( PairingAlgo::Mortar ),
          _cont_init( InitialState::Interpenetrated ),
          _dist_ratio( -1.0 ),

          _beam( false ),
          _dist_supp( nullptr ),
          _shell( false ),
          _cara( nullptr ) {};

    /** @brief restricted constructor (Set) and method (Get) to support pickling */
    PairingParameter( const py::tuple &tup );
    py::tuple _getState() const;

    PairingAlgo getAlgorithm() const { return _algo; };

    ASTERDOUBLE getDistanceRatio() const { return _dist_ratio; };

    InitialState getInitialState() const { return _cont_init; };

    GenericFunctionPtr getDistanceFunction() const { return _dist_supp; };

    ElementaryCharacteristicsPtr getElementaryCharacteristics() const { return _cara; };

    void setAlgorithm( const PairingAlgo &algo ) { _algo = algo; };

    void setDistanceRatio( const ASTERDOUBLE &dist_ratio ) { _dist_ratio = dist_ratio; };

    void setInitialState( const InitialState &cont_init ) { _cont_init = cont_init; };

    void setDistanceFunction( const GenericFunctionPtr &dist_supp ) { _dist_supp = dist_supp; };

    void enableBeamDistance( const bool &beam ) { _beam = beam; }

    bool hasBeamDistance() const { return _beam; }

    void enableShellDistance( const bool &shell ) { _shell = shell; }

    bool hasShellDistance() const { return _shell; }

    void setElementaryCharacteristics( const ElementaryCharacteristicsPtr &cara ) { _cara = cara; };
};

/**
 * @typedef PairingParameterPtr
 * @brief Pointeur intelligent vers un PairingParameter
 */
typedef std::shared_ptr< PairingParameter > PairingParameterPtr;

#endif /* CONTACT_PARAM_H_ */
