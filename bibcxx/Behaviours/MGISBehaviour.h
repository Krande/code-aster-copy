/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

#pragma once

#include "astercxx.h"

#include "DataStructures/DataStructure.h"
#include "Supervis/ResultNaming.h"

#include <iostream>
#include <sstream>

#ifdef ASTER_HAVE_MGIS

#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"

namespace MGB = mgis::behaviour;

using BehaviourPtr = std::shared_ptr< MGB::Behaviour >;
using BehaviourDataPtr = std::shared_ptr< MGB::BehaviourData >;

enum struct StateId { s0, s1 };

class MGISBehaviour : public DataStructure {
  private:
    BehaviourPtr _behav;
    BehaviourDataPtr _data;
    std::string _libpath;
    std::string _bname;
    VectorReal _rotmatrix, _rotmatrixT;
    bool _initialized;
    JeveuxVectorChar16 _addr;

  public:
    using MGISBehaviourPtr = std::shared_ptr< MGISBehaviour >;

    MGISBehaviour( const std::string name )
        : DataStructure( name, 8, "COMPOR_MGIS" ),
          _behav( nullptr ),
          _data( nullptr ),
          _libpath( "" ),
          _bname( "" ),
          _initialized( false ),
          _addr( JeveuxVectorChar16( getName() + ".ADDR" ) ) {
        if ( !_addr->exists() ) {
            _addr->allocate( 1 );
        } else {
            _addr->updateValuePointer();
        }
        std::stringstream sstr;
        sstr << std::hex << id();
        ( *_addr )[0] = sstr.str();
#ifdef ASTER_DEBUG_CXX
        std::cout << "creating MGISBehaviour with id: " << id() << " (0x" << sstr.str() << ")"
                  << std::endl;
#endif
    };

    MGISBehaviour() : MGISBehaviour( ResultNaming::getNewResultName() ) {};

    ~MGISBehaviour() { std::cout << "deleting MGISBehaviour with id: " << id() << std::endl; }

    void setLibPath( std::string path ) { _libpath = path; };
    void setBehaviourName( std::string behav_name ) { _bname = behav_name; };
    VectorString internal_state() const;

    MGB::State *getState( const StateId stid ) const;
    void repr() const;

    /* External State VariableS */
    int getNumberOfExternalStateVariables() const;
    VectorString getExternalStateVariablesNames() const;
    void setExternalStateVariables( const StateId stid, const ASTERDOUBLE *values );

    /* Internal State VariableS */
    int getNumberOfInternalStateVariables() const;
    int getSizeOfInternalStateVariables() const;
    VectorString getInternalStateVariablesNames() const;
    VectorLong getInternalStateVariablesTypes() const;
    VectorLong getInternalStateVariablesSizes() const;
    void setInternalStateVariables( const StateId stid, const ASTERDOUBLE *values );

    /* Other State Variables */
    int getSizeOfGradients() const;
    void setGradients( const StateId stid, const ASTERDOUBLE *values, const int insize );
    int getSizeOfThermodynamicForces() const;
    void setThermodynamicForces( const StateId stid, const ASTERDOUBLE *values, const int insize );
    void setRotationMatrix( const ASTERDOUBLE *drot );

    /* Parameters */
    void setParameter( const std::string param, const ASTERDOUBLE value );
    void setParameter( const std::string param, const ASTERINTEGER value );

    /// @brief Get value of keyword from MFront
    /// @tparam T Type of return (ASTERDOUBLE or ASTERINTEGER)
    /// @tparam U Type for cast (double for ASTERDOUBLE, int for ASTERINTEGER)
    /// @param param keyword parameter
    /// @return value
    template < typename T, typename U >
    T getMFrontParameter( const std::string param );

    /* Material Properties */
    int getNumberOfMaterialProperties() const;
    VectorString getMaterialPropertiesNames() const;
    void setMaterialProperties( const StateId stid, const ASTERDOUBLE *values );

    /* General */
    void _load_library( const MGB::Hypothesis hyp, const bool finiteStrain = false,
                        const MGB::FiniteStrainBehaviourOptions::StressMeasure stress_measure =
                            MGB::FiniteStrainBehaviourOptions::CAUCHY,
                        const MGB::FiniteStrainBehaviourOptions::TangentOperator tangent_op =
                            MGB::FiniteStrainBehaviourOptions::DSIG_DF );
    void load( const MGB::Hypothesis hyp );
    void load( const MGB::Hypothesis hyp,
               const MGB::FiniteStrainBehaviourOptions::StressMeasure stress_measure,
               const MGB::FiniteStrainBehaviourOptions::TangentOperator tangent_op );
    void setInitialState( ASTERDOUBLE rdt, ASTERDOUBLE dt, ASTERDOUBLE *K );
    int integrate();
    void setOutputs( ASTERDOUBLE *stress, ASTERDOUBLE *internal_variables, ASTERDOUBLE *K,
                     ASTERDOUBLE *rdt );
};

void MGISExceptionHandler();

using MGISBehaviourPtr = std::shared_ptr< MGISBehaviour >;

void convertTensorFromMgis( const VectorReal &tens_vec, ASTERDOUBLE *array );

void convertMatrixFromMgis( const VectorReal &tens_vec, ASTERDOUBLE *array );

void convertTensorToMgis( const ASTERDOUBLE *src, const int insize, VectorReal &dest,
                          const int destsize );

#endif // ASTER_HAVE_MGIS
