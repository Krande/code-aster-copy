#ifndef LINESEARCHMETHOD_H_
#define LINESEARCHMETHOD_H_

/**
 * @file LineSearchMethod.h
 * @brief Definition of the linesearch method
 * @author Natacha Béreux
 * @section LICENCE
 *   Copyright (C) 1991 - 2020  EDF R&D                www.code-aster.org
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

/* person_in_charge: natacha.bereux at edf.fr */
#include "astercxx.h"

#include "Solvers/SolverControl.h"
#include "Utilities/GenericParameter.h"

enum LineSearchEnum { Corde, Mixte, Pilotage };
const int nbLineSearch = 4;
extern const char *LineSearchNames[nbLineSearch];

class LineSearchMethodClass {
  private:
    /** @brief LineSearch Method */
    LineSearchEnum _lineSearchMethod;
    /** @brief Contrôle de la convergence de la méthode  */
    SolverControlPtr _control;
    /** LineSearch method name */
    GenParamPtr _methode;
    /** Intervalle de recherche de rho */
    GenParamPtr _rhoMin;
    GenParamPtr _rhoMax;
    GenParamPtr _rhoExcl;
    /** Control */
    GenParamPtr _resi_line_rela;
    GenParamPtr _iter_line_maxi;

    ListGenParam _listOfParameters;

  public:
    /**
     * @brief Constructeur
     */
    LineSearchMethodClass( LineSearchEnum curLineSearch = Corde )
        : _lineSearchMethod( curLineSearch ),
          _control( SolverControlPtr( new SolverControlClass() ) ),
          _methode(boost::make_shared<GenParam>( "METHODE", false )),
          _rhoMin(boost::make_shared<GenParam>( "RHO_MIN", 1.e-2, false )),
          _rhoMax(boost::make_shared<GenParam>( "RHO_MAX", 1.e1, false )),
          _rhoExcl(boost::make_shared<GenParam>( "RHO_EXCL", 9.e-3, false )),
          _resi_line_rela(boost::make_shared<GenParam>( "RESI_LINE_RELA", false )),
          _iter_line_maxi(boost::make_shared<GenParam>( "ITER_LINE_MAXI", false )) {
        _control->setRelativeTolerance( 1.e-1 );
        _control->setMaximumNumberOfIterations( 3 );

        _methode->setValue( std::string( LineSearchNames[(int)curLineSearch]) );

        _resi_line_rela->setValue( _control->getRelativeTolerance());
        _iter_line_maxi->setValue( _control->getMaximumNumberOfIterations());

        _listOfParameters.push_back( _methode );
        _listOfParameters.push_back( _rhoMin );
        _listOfParameters.push_back( _rhoMax );
        _listOfParameters.push_back( _rhoExcl );
        _listOfParameters.push_back( _resi_line_rela );
        _listOfParameters.push_back( _iter_line_maxi );
    };
    /**
    @brief set minimum value of rho
    */
    void setMinimumRhoValue( double rhoMin ) { _rhoMin->setValue( rhoMin); };
    /**
    @brief set maximum value of rho
    */
    void setMaximumRhoValue( double rhoMax ) { _rhoMax->setValue( rhoMax); };
    /**
    @brief
    */
    void setExclRhoValue( double rhoExcl ) { _rhoExcl->setValue( rhoExcl); };
    /**
    @brief set maximum number of iterations
    */
    void setMaximumNumberOfIterations( ASTERINTEGER nIterMax ) {
        _iter_line_maxi->setValue( nIterMax);
        _control->setMaximumNumberOfIterations( nIterMax );
    };
    /**
    @brief set the relative tolerance
    */
    void setRelativeTolerance( double reslin ) {
        _resi_line_rela->setValue( reslin);
        _control->setRelativeTolerance( reslin );
    };
    /**
     * @brief Récupération de la liste des paramètres
     * @return Liste constante des paramètres déclarés
     */
    const ListGenParam &getListOfParameters() { return _listOfParameters; };
};

/**
 * @typedef NonLinearMethodPtr
 * @brief Enveloppe d'un pointeur intelligent vers un LineSearchMethodClass
 */
typedef boost::shared_ptr< LineSearchMethodClass > LineSearchMethodPtr;

#endif /* LINESEARCHMETHOD_H_ */
