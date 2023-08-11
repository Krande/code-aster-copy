/**
 * @file DiscreteComputation.cxx
 * @brief Implementation of class DiscreteComputation
 * @section LICENCE
 *   Copyright (C) 1991 2023  EDF R&D                www.code-aster.org
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

#include "aster_fort_calcul.h"
#include "aster_fort_superv.h"

#include "Discretization/Calcul.h"
#include "Discretization/DiscreteComputation.h"
#include "Loads/DirichletBC.h"
#include "Loads/MechanicalLoad.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Modeling/XfemModel.h"
#include "Utilities/Tools.h"

std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getThermalNeumannForces( const ASTERDOUBLE time_curr,
                                              const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    auto elemVect = std::make_shared< ElementaryVectorTemperatureReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_THER" );
    elemVect->prepareCompute( calcul_option );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto listOfLoads = _phys_problem->getListOfLoads();
    auto model_FEDesc = currModel->getFiniteElementDescriptor();
    AS_ASSERT( model_FEDesc );
    auto isXfem = currModel->existsXfem();

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );

    auto therLoadReal = listOfLoads->getThermalLoadsReal();
    for ( const auto &load : therLoadReal ) {
        // Termes FLUX XYZ
        if ( load->hasLoadField( "FLURE" ) ) {
            auto flow_xyz_field = load->getConstantLoadField( "FLURE" );
            calcul->setOption( "CHAR_THER_FLUN_R" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addTimeField( "PTEMPSR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PFLUXNR", flow_xyz_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        // Termes FLUX NORM
        if ( load->hasLoadField( "FLUR2" ) ) {
            auto flow_nor_field = load->getConstantLoadField( "FLUR2" );
            calcul->setOption( "CHAR_THER_FLUX_R" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addTimeField( "PTEMPSR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PFLUXVR", flow_nor_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        iload++;
    }

    auto therLoadFunc = listOfLoads->getThermalLoadsFunction();
    for ( const auto &load : therLoadFunc ) {
        // Termes FLUX XYZ
        if ( load->hasLoadField( "FLURE" ) ) {
            auto flow_xyz_field = load->getConstantLoadField( "FLURE" );
            calcul->setOption( "CHAR_THER_FLUN_F" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addTimeField( "PTEMPSR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PFLUXNF", flow_xyz_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        // Termes FLUX NORM
        if ( load->hasLoadField( "FLUR2" ) ) {
            auto flow_nor_field = load->getConstantLoadField( "FLUR2" );
            calcul->setOption( "CHAR_THER_FLUX_F" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addTimeField( "PTEMPSR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PFLUXVF", flow_nor_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        iload++;
    }

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            return elemVect->assembleWithLoadFunctions( _phys_problem->getDOFNumbering(),
                                                        time_curr );
        } else {
            FieldOnNodesRealPtr vectAsse =
                std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
};

std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getThermalVolumetricForces( const ASTERDOUBLE time_curr,
                                                 const FieldOnCellsRealPtr varc_curr,
                                                 const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    auto elemVect = std::make_shared< ElementaryVectorTemperatureReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_THER" );
    elemVect->prepareCompute( calcul_option );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto listOfLoads = _phys_problem->getListOfLoads();
    auto model_FEDesc = currModel->getFiniteElementDescriptor();
    AS_ASSERT( model_FEDesc );
    auto isXfem = currModel->existsXfem();

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );

    auto therLoadReal = listOfLoads->getThermalLoadsReal();
    for ( const auto &load : therLoadReal ) {

        // Termes SOURCE
        if ( load->hasLoadField( "SOURE" ) ) {
            auto source_field = load->getConstantLoadField( "SOURE" );
            calcul->setOption( "CHAR_THER_SOUR_R" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addTimeField( "PTEMPSR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PSOURCR", source_field );

            if ( currMater && currMater->hasExternalStateVariable() ) {
                if ( !varc_curr || !varc_curr->exists() ) {
                    raiseAsterError( "External state variables are needed but not given" );
                }
                calcul->addInputField( "PVARCPR", varc_curr );
            }

            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }
        // Termes SOURCE CALCULEE
        if ( load->hasLoadField( "SOURC" ) ) {
            auto computed_source_field = load->getLoadField( "SOURC" );
            calcul->setOption( "CHAR_THER_SOUR_R" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addTimeField( "PTEMPSR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PSOURCR", computed_source_field );

            if ( currMater && currMater->hasExternalStateVariable() ) {
                if ( !varc_curr || !varc_curr->exists() ) {
                    raiseAsterError( "External state variables are needed but not given" );
                }
                calcul->addInputField( "PVARCPR", varc_curr );
            }

            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        // Termes PRE_GRAD_TEMP
        if ( load->hasLoadField( "GRAIN" ) ) {
            auto pregrad_field = load->getConstantLoadField( "GRAIN" );
            calcul->setOption( "CHAR_THER_GRAI_R" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            if ( currMater ) {
                calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

                if ( currMater->hasExternalStateVariable() ) {
                    if ( !varc_curr || !varc_curr->exists() ) {
                        raiseAsterError( "External state variables are needed but not given" );
                    }
                    calcul->addInputField( "PVARCPR", varc_curr );
                }
            }
            calcul->addInputField( "PGRAINR", pregrad_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        iload++;
    }

    auto therLoadFunc = listOfLoads->getThermalLoadsFunction();
    for ( const auto &load : therLoadFunc ) {
        auto load_FEDesc = load->getFiniteElementDescriptor();

        if ( load->hasLoadField( "SOURE" ) ) {
            auto source_field = load->getConstantLoadField( "SOURE" );
            calcul->setOption( "CHAR_THER_SOUR_F" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addTimeField( "PTEMPSR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PSOURCF", source_field );

            if ( currMater && currMater->hasExternalStateVariable() ) {
                if ( !varc_curr || !varc_curr->exists() ) {
                    raiseAsterError( "External state variables are needed but not given" );
                }
                calcul->addInputField( "PVARCPR", varc_curr );
            }

            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        // Termes PRE_GRAD_TEMP
        if ( load->hasLoadField( "GRAIN" ) ) {
            auto pregrad_field = load->getConstantLoadField( "GRAIN" );
            calcul->setOption( "CHAR_THER_GRAI_F" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addTimeField( "PTEMPSR", time_curr, 0.0, -1.0 );

            if ( currMater ) {
                calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

                if ( currMater->hasExternalStateVariable() ) {
                    if ( !varc_curr || !varc_curr->exists() ) {
                        raiseAsterError( "External state variables are needed but not given" );
                    }
                    calcul->addInputField( "PVARCPR", varc_curr );
                }
            }
            calcul->addInputField( "PGRAINF", pregrad_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        iload++;
    }

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            return elemVect->assembleWithLoadFunctions( _phys_problem->getDOFNumbering(),
                                                        time_curr );
        } else {
            FieldOnNodesRealPtr vectAsse =
                std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
};

std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getThermalExchangeForces( const FieldOnNodesRealPtr temp_curr,
                                               const ASTERDOUBLE time_curr,
                                               const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    auto elemVect = std::make_shared< ElementaryVectorTemperatureReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_THER" );
    elemVect->prepareCompute( calcul_option );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto listOfLoads = _phys_problem->getListOfLoads();
    auto model_FEDesc = currModel->getFiniteElementDescriptor();
    AS_ASSERT( model_FEDesc );
    auto isXfem = currModel->existsXfem();

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );

    auto therLoadReal = listOfLoads->getThermalLoadsReal();
    for ( const auto &load : therLoadReal ) {
        auto load_FEDesc = load->getFiniteElementDescriptor();

        // Termes ECHANGE
        if ( load->hasLoadField( "COEFH" ) && load->hasLoadField( "T_EXT" ) ) {
            auto exchange_field = load->getConstantLoadField( "COEFH" );
            auto ext_temp_field = load->getConstantLoadField( "T_EXT" );

            calcul->setOption( "CHAR_THER_TEXT_R" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addInputField( "PTEMPER", temp_curr );
            calcul->addTimeField( "PTEMPSR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PCOEFHR", exchange_field );
            calcul->addInputField( "PT_EXTR", ext_temp_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        // Termes ECHANGE_PAROI
        if ( load->hasLoadField( "HECHP" ) ) {
            auto wall_exchange_field = load->getConstantLoadField( "HECHP" );
            calcul->setOption( "CHAR_THER_PARO_R" );
            calcul->clearInputs();
            if ( isXfem ) {
                XfemModelPtr currXfemModel = currModel->getXfemModel();
                calcul->addXFEMField( currXfemModel );
                calcul->setFiniteElementDescriptor( model_FEDesc );
            } else {
                calcul->setFiniteElementDescriptor( load_FEDesc );
            }
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addInputField( "PTEMPER", temp_curr );
            calcul->addTimeField( "PTEMPSR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PHECHPR", wall_exchange_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        iload++;
    }

    auto therLoadFunc = listOfLoads->getThermalLoadsFunction();
    for ( const auto &load : therLoadFunc ) {
        auto load_FEDesc = load->getFiniteElementDescriptor();
        // Termes ECHANGE
        if ( load->hasLoadField( "COEFH" ) && load->hasLoadField( "T_EXT" ) ) {
            auto exchange_field = load->getConstantLoadField( "COEFH" );
            auto ext_temp_field = load->getConstantLoadField( "T_EXT" );

            calcul->setOption( "CHAR_THER_TEXT_F" );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addInputField( "PTEMPER", temp_curr );
            calcul->addTimeField( "PTEMPSR", time_curr, 0.0, -1.0 );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->addInputField( "PCOEFHF", exchange_field );
            calcul->addInputField( "PT_EXTF", ext_temp_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        // Termes ECHANGE_PAROI
        if ( load->hasLoadField( "HECHP" ) ) {
            auto wall_exchange_field = load->getConstantLoadField( "HECHP" );
            calcul->setOption( "CHAR_THER_PARO_F" );
            calcul->clearInputs();
            if ( isXfem ) {
                XfemModelPtr currXfemModel = currModel->getXfemModel();
                calcul->addXFEMField( currXfemModel );
                calcul->setFiniteElementDescriptor( model_FEDesc );
            } else {
                calcul->setFiniteElementDescriptor( load_FEDesc );
            }
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addInputField( "PTEMPER", temp_curr );
            calcul->addTimeField( "PTEMPSR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PHECHPF", wall_exchange_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }
        iload++;
    }

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            return elemVect->assembleWithLoadFunctions( _phys_problem->getDOFNumbering(),
                                                        time_curr );
        } else {
            FieldOnNodesRealPtr vectAsse =
                std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
};

std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getThermalNonLinearNeumannForces( const FieldOnNodesRealPtr temp_curr,
                                                       const ASTERDOUBLE time_curr,
                                                       const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    auto elemVect = std::make_shared< ElementaryVectorTemperatureReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_THER" );
    elemVect->prepareCompute( calcul_option );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto listOfLoads = _phys_problem->getListOfLoads();
    auto model_FEDesc = currModel->getFiniteElementDescriptor();
    AS_ASSERT( model_FEDesc );

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );

    auto impl = [&]( auto load, const ASTERINTEGER &load_i, const std::string &option,
                     const std::string &name, const std::string &param,
                     const FiniteElementDescriptorPtr FED ) {
        if ( load->hasLoadField( name ) ) {
            calcul->setOption( option );
            calcul->setFiniteElementDescriptor( FED );

            calcul->clearInputs();
            calcul->addTimeField( "PTEMPSR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addInputField( "PTEMPER", temp_curr );

            calcul->addInputField( param, load->getConstantLoadField( name ) );

            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             load_i );
            }
        }
    };

    auto therLoadReal = listOfLoads->getThermalLoadsReal();
    for ( const auto &load : therLoadReal ) {
        impl( load, iload, "CHAR_THER_RAYO_R", "RAYO", "PRAYONR", model_FEDesc );

        iload++;
    }

    auto therLoadFunc = listOfLoads->getThermalLoadsFunction();
    for ( const auto &load : therLoadFunc ) {
        impl( load, iload, "CHAR_THER_FLUNL", "FLUNL", "PFLUXNL", model_FEDesc );
        impl( load, iload, "CHAR_THER_RAYO_F", "RAYO", "PRAYONF", model_FEDesc );

        iload++;
    }

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            return elemVect->assembleWithLoadFunctions( _phys_problem->getDOFNumbering(),
                                                        time_curr );
        } else {
            FieldOnNodesRealPtr vectAsse =
                std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
};

std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getThermalNonLinearVolumetricForces( const FieldOnNodesRealPtr temp_curr,
                                                          const ASTERDOUBLE time_curr,
                                                          const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    auto elemVect = std::make_shared< ElementaryVectorTemperatureReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_THER" );
    elemVect->prepareCompute( calcul_option );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto listOfLoads = _phys_problem->getListOfLoads();
    auto model_FEDesc = currModel->getFiniteElementDescriptor();
    AS_ASSERT( model_FEDesc );

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );

    auto impl = [&]( auto load, const ASTERINTEGER &load_i, const std::string &option,
                     const std::string &name, const std::string &param,
                     const FiniteElementDescriptorPtr FED ) {
        if ( load->hasLoadField( name ) ) {
            calcul->setOption( option );
            calcul->setFiniteElementDescriptor( FED );

            calcul->clearInputs();
            calcul->addTimeField( "PTEMPSR", time_curr, 0.0, 0.0 );
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addInputField( "PTEMPER", temp_curr );

            calcul->addInputField( param, load->getConstantLoadField( name ) );

            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             load_i );
            }
        }
    };

    auto therLoadReal = listOfLoads->getThermalLoadsReal();
    for ( const auto &load : therLoadReal ) {
        impl( load, iload, "CHAR_THER_SOURNL", "SOUNL", "PSOURNL", model_FEDesc );

        iload++;
    }

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            return elemVect->assembleWithLoadFunctions( _phys_problem->getDOFNumbering(),
                                                        time_curr );
        } else {
            FieldOnNodesRealPtr vectAsse =
                std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
};

std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getTransientThermalLoadForces( const ASTERDOUBLE time_curr,
                                                    const ASTERDOUBLE time_step,
                                                    const ASTERDOUBLE theta,
                                                    const FieldOnNodesRealPtr _previousPrimalField,
                                                    const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    auto elemVect = std::make_shared< ElementaryVectorTemperatureReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_THER" );
    elemVect->prepareCompute( calcul_option );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto listOfLoads = _phys_problem->getListOfLoads();
    auto model_FEDesc = currModel->getFiniteElementDescriptor();
    AS_ASSERT( model_FEDesc );

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );

    auto therLoadReal = listOfLoads->getThermalLoadsReal();
    for ( const auto &load : therLoadReal ) {
        if ( load->hasLoadResult() ) {
            std::string evol_char_name = load->getLoadResultName();
            std::string para_flun( "FLUN" );
            std::string para_coefh( "COEF_H" );
            std::string para_text( "T_EXT" );
            std::string access_var( "INST" );
            std::string base( "G" );
            std::string extr_right( "EXCLU" );
            std::string extr_left( "EXCLU" );
            ASTERINTEGER iret = 100;
            ASTERINTEGER stop = 0;

            FieldOnCellsRealPtr evol_flow_xyz_field =
                std::make_shared< FieldOnCellsReal >( model_FEDesc );
            // On cherche le champ FLUN. Si il existe on calcule l'option CHAR_THER_FLUN_R
            // Si il n'existe pas on suppose l'existence des champs pour calculer CHAR_THER_TEXT_R
            CALLO_RSINCH( evol_char_name, para_flun, access_var, &time_curr,
                          evol_flow_xyz_field->getName(), extr_right, extr_left, &stop, base,
                          &iret );

            if ( iret >= 2 ) {

                FieldOnCellsRealPtr evol_exchange_field =
                    std::make_shared< FieldOnCellsReal >( model_FEDesc );
                FieldOnCellsRealPtr evol_ext_temp_field =
                    std::make_shared< FieldOnCellsReal >( model_FEDesc );

                CALLO_RSINCH( evol_char_name, para_coefh, access_var, &time_curr,
                              evol_exchange_field->getName(), extr_right, extr_left, &stop, base,
                              &iret );

                if ( iret >= 2 ) {
                    AS_ABORT( "Cannot find COEF_H in EVOL_CHAR " + evol_char_name + " at time " +
                              std::to_string( time_curr ) );
                }

                CALLO_RSINCH( evol_char_name, para_text, access_var, &time_curr,
                              evol_ext_temp_field->getName(), extr_right, extr_left, &stop, base,
                              &iret );

                if ( iret >= 2 ) {
                    AS_ABORT( "Cannot find T_EXT in EVOL_CHAR " + evol_char_name + " at time " +
                              std::to_string( time_curr ) );
                }
                AS_ASSERT( _previousPrimalField && _previousPrimalField->exists() );

                calcul->setOption( "CHAR_THER_TEXT_R" );
                calcul->setFiniteElementDescriptor( model_FEDesc );
                calcul->clearInputs();
                calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
                calcul->addInputField( "PTEMPER", _previousPrimalField );
                calcul->addTimeField( "PTEMPSR", time_curr, time_step, theta );
                calcul->addInputField( "PCOEFHR", evol_exchange_field );
                calcul->addInputField( "PT_EXTR", evol_ext_temp_field );
                calcul->clearOutputs();
                calcul->addOutputElementaryTerm( "PVECTTR",
                                                 std::make_shared< ElementaryTermReal >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                    elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                                 iload );
                }
            } else {
                calcul->setOption( "CHAR_THER_FLUN_R" );
                calcul->setFiniteElementDescriptor( model_FEDesc );
                calcul->clearInputs();
                calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
                calcul->addTimeField( "PTEMPSR", time_curr, time_step, theta );
                calcul->addInputField( "PFLUXNR", evol_flow_xyz_field );
                calcul->clearOutputs();
                calcul->addOutputElementaryTerm( "PVECTTR",
                                                 std::make_shared< ElementaryTermReal >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                    elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                                 iload );
                }
            }
        }

        iload++;
    }

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            return elemVect->assembleWithLoadFunctions( _phys_problem->getDOFNumbering(),
                                                        time_curr );
        } else {
            FieldOnNodesRealPtr vectAsse =
                std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
};

/** @brief Compute AFFE_CHAR_THER TEMP_IMPO */
std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getThermalImposedDualBC( const ASTERDOUBLE time_curr,
                                              const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    auto elemVect = std::make_shared< ElementaryVectorTemperatureReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_THER" );
    elemVect->prepareCompute( calcul_option );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto listOfLoads = _phys_problem->getListOfLoads();

    auto calcul = std::make_unique< Calcul >( calcul_option );

    auto impl = [&]( auto loads, bool real ) {
        std::string name;
        if ( real ) {
            calcul->setOption( "THER_DDLI_R" );
            name = "PDDLIMR";
        } else {
            calcul->setOption( "THER_DDLI_F" );
            name = "PDDLIMF";
        }
        for ( const auto &load : loads ) {
            auto load_FEDesc = load->getFiniteElementDescriptor();
            auto impo_field = load->getImposedField();
            if ( impo_field && impo_field->exists() && load_FEDesc ) {
                calcul->clearInputs();
                calcul->clearOutputs();
                calcul->setFiniteElementDescriptor( load_FEDesc );
                calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
                calcul->addInputField( name, impo_field );
                calcul->addTimeField( "PTEMPSR", time_curr, 0.0, -1.0 );
                calcul->addOutputElementaryTerm( "PVECTTR",
                                                 std::make_shared< ElementaryTermReal >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                    elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                                 iload );
                }
            }
            iload++;
        }
    };

    impl( listOfLoads->getThermalLoadsReal(), true );
    impl( listOfLoads->getThermalLoadsFunction(), false );

#ifdef ASTER_HAVE_MPI
    impl( listOfLoads->getParallelThermalLoadsReal(), true );
    impl( listOfLoads->getParallelThermalLoadsFunction(), false );
#endif

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            return elemVect->assembleWithLoadFunctions( _phys_problem->getDOFNumbering(),
                                                        time_curr );
        } else {
            FieldOnNodesRealPtr vectAsse =
                std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
}

/** @brief Compute CHAR_THER_EVOL */
FieldOnNodesRealPtr DiscreteComputation::getTransientThermalForces(
    const ASTERDOUBLE time_curr, const ASTERDOUBLE time_step, const ASTERDOUBLE theta,
    const FieldOnNodesRealPtr previousPrimalField, const FieldOnCellsRealPtr varc_curr ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );
    AS_ASSERT( previousPrimalField && previousPrimalField->exists() );

    auto elemVect = std::make_shared< ElementaryVectorReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );

    // Setup
    const std::string calcul_option( "CHAR_THER_EVOL" );
    elemVect->prepareCompute( calcul_option );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );
    calcul->clearInputs();
    calcul->clearOutputs();

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    calcul->addTimeField( "PTEMPSR", time_curr, time_step, theta );
    calcul->addInputField( "PTEMPER", previousPrimalField );

    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

        if ( currMater->hasExternalStateVariable() ) {
            if ( !varc_curr || !varc_curr->exists() ) {
                raiseAsterError( "External state variables are needed but not given" );
            }
            calcul->addInputField( "PVARCPR", varc_curr );
        }
    }

    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
    }
    if ( currModel->existsXfem() ) {
        XfemModelPtr currXfemModel = currModel->getXfemModel();
        calcul->addXFEMField( currXfemModel );
    }
    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );

    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) )
            elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ) );
    };

    elemVect->build();
    // Assemble
    return elemVect->assemble( _phys_problem->getDOFNumbering() );
}

/** @brief Compute CHAR_THER_EVOLNI */
FieldOnNodesRealPtr DiscreteComputation::getNonLinearTransientThermalForces(
    const FieldOnNodesRealPtr temp_prev, const FieldOnNodesRealPtr temp_step,
    const ASTERDOUBLE time_prev, const ASTERDOUBLE time_step, const ASTERDOUBLE theta,
    const FieldOnCellsRealPtr varc_curr ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    auto elemVect = std::make_shared< ElementaryVectorReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );

    // Setup
    const std::string calcul_option( "CHAR_THER_EVOLNI" );
    elemVect->prepareCompute( calcul_option );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currBehav = _phys_problem->getBehaviourProperty();

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );
    calcul->clearInputs();
    calcul->clearOutputs();

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    calcul->addTimeField( "PTEMPSR", time_prev + time_step, time_step, theta );
    auto temp_curr = std::make_shared< FieldOnNodesReal >( *temp_prev + *temp_step );
    calcul->addInputField( "PTEMPER", temp_prev );

    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

        if ( currMater->hasExternalStateVariable() ) {
            if ( !varc_curr || !varc_curr->exists() ) {
                raiseAsterError( "External state variables are needed but not given" );
            }
            calcul->addInputField( "PVARCPR", varc_curr );
        }
    }

    if ( currBehav ) {
        calcul->addInputField( "PCOMPOR", currBehav->getBehaviourField() );
    }

    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
    }
    if ( currModel->existsXfem() ) {
        XfemModelPtr currXfemModel = currModel->getXfemModel();
        calcul->addXFEMField( currXfemModel );
    }
    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PVECTTI", std::make_shared< ElementaryTermReal >() );
    calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );

    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PVECTTI" ) )
            elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTI" ) );
        // Do not add linear part
        // if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) )
        //     elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ) );
    };
    elemVect->build();
    // Assemble
    return elemVect->assemble( _phys_problem->getDOFNumbering() );
}

/**
 * @brief Compute elementary forces for internal forces (RAPH_THER)
 */
FieldOnNodesRealPtr
DiscreteComputation::getInternalThermalForces( const FieldOnNodesRealPtr temp_step,
                                               const FieldOnCellsRealPtr varc_curr,
                                               const VectorString &groupOfCells ) const {
    AS_ASSERT( _phys_problem->getModel()->isThermal() );
    const std::string option( "RAPH_THER" );

    auto elemVect = std::make_shared< ElementaryVectorReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );
    elemVect->prepareCompute( option );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currBehav = _phys_problem->getBehaviourProperty();

    // Prepare computing
    auto calcul = std::make_shared< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );

    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

        if ( currMater->hasExternalStateVariable() ) {
            if ( !varc_curr || !varc_curr->exists() ) {
                raiseAsterError( "External state variables are needed but not given" );
            }
            calcul->addInputField( "PVARCPR", varc_curr );
        }
    }

    if ( currBehav ) {
        calcul->addInputField( "PCOMPOR", currBehav->getBehaviourField() );
    }

    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
    }

    if ( currModel->existsXfem() ) {
        XfemModelPtr currXfemModel = currModel->getXfemModel();
        calcul->addXFEMField( currXfemModel );
    }

    // Add Thermal Field
    calcul->addInputField( "PTEMPEI", temp_step );

    // TODO:
    // calcul->addInputField( "PTMPCHF", dry_curr );

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PRESIDU", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for mass
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PRESIDU" ) )
            elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PRESIDU" ) );
    };

    elemVect->build();
    return elemVect->assemble( _phys_problem->getDOFNumbering() );
}

/**
 * @brief Compute elementary forces for capacitiy forces (MASS_THER_RESI)
 */
FieldOnNodesRealPtr DiscreteComputation::getNonLinearCapacityForces(
    const FieldOnNodesRealPtr temp_prev, const FieldOnNodesRealPtr temp_step,
    const ASTERDOUBLE &time_step, const FieldOnCellsRealPtr varc_curr,
    const VectorString &groupOfCells ) const {
    AS_ASSERT( _phys_problem->getModel()->isThermal() );
    const std::string option( "MASS_THER_RESI" );

    auto elemVect = std::make_shared< ElementaryVectorReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );
    elemVect->prepareCompute( option );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currBehav = _phys_problem->getBehaviourProperty();

    // Prepare computing
    auto calcul = std::make_shared< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );

    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

        if ( currMater->hasExternalStateVariable() ) {
            if ( !varc_curr || !varc_curr->exists() ) {
                raiseAsterError( "External state variables are needed but not given" );
            }
            calcul->addInputField( "PVARCPR", varc_curr );
        }
    }

    if ( currBehav ) {
        calcul->addInputField( "PCOMPOR", currBehav->getBehaviourField() );
    }

    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
    }

    if ( currModel->existsXfem() ) {
        XfemModelPtr currXfemModel = currModel->getXfemModel();
        calcul->addXFEMField( currXfemModel );
    }

    // Add time field
    calcul->addTimeField( "PTEMPSR", 0.0, time_step, 0.0 );

    // Current Thermal Field
    calcul->addInputField( "PTEMPER", temp_prev );
    calcul->addInputField( "PTEMPEI", temp_step );

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PRESIDU", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for mass
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PRESIDU" ) )
            elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PRESIDU" ) );
    };

    elemVect->build();
    return elemVect->assemble( _phys_problem->getDOFNumbering() );
}

FieldOnNodesRealPtr DiscreteComputation::getDualTemperature( FieldOnNodesRealPtr temp_curr,
                                                             ASTERDOUBLE scaling ) const {

    auto elemVect = std::make_shared< ElementaryVectorReal >(
        _phys_problem->getModel(), _phys_problem->getMaterialField(),
        _phys_problem->getElementaryCharacteristics(), _phys_problem->getListOfLoads() );

    elemVect->prepareCompute( "THER_BU_R" );

    if ( !_phys_problem->isThermal() ) {
        AS_ABORT( "Not implemented" );
    };

    // Prepare loads
    const auto listOfLoads = _phys_problem->getListOfLoads();

    auto calcul = std::make_unique< Calcul >( "THER_BU_R" );

    // ConstantField for scaling
    auto cf_scaling = std::make_shared< ConstantFieldOnCellsReal >( _phys_problem->getMesh() );
    cf_scaling->allocate( "NEUT_R" );
    ConstantFieldOnZone a( _phys_problem->getMesh() );
    ConstantFieldValues< ASTERDOUBLE > b( {"X1"}, {scaling} );
    cf_scaling->setValueOnZone( a, b );

    auto impl = [&]( auto loads ) {
        for ( const auto &load : loads ) {
            auto FEDesc = load->getFiniteElementDescriptor();
            auto field = load->getMultiplicativeField();
            if ( field && field->exists() && FEDesc ) {
                calcul->clearInputs();
                calcul->clearOutputs();
                calcul->setFiniteElementDescriptor( FEDesc );
                calcul->addInputField( "PDDLMUR", field );
                calcul->addInputField( "PDDLIMR", temp_curr );
                calcul->addInputField( "PALPHAR", cf_scaling );
                calcul->addOutputElementaryTerm( "PVECTTR",
                                                 std::make_shared< ElementaryTermReal >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                    elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ) );
                }
            }
        }
    };

    impl( listOfLoads->getThermalLoadsReal() );
    impl( listOfLoads->getThermalLoadsFunction() );

#ifdef ASTER_HAVE_MPI
    impl( listOfLoads->getParallelThermalLoadsReal() );
    impl( listOfLoads->getParallelThermalLoadsFunction() );
#endif

    // Assemble
    elemVect->build();

    FieldOnNodesRealPtr buth;
    if ( elemVect->hasElementaryTerm() ) {
        buth = elemVect->assemble( _phys_problem->getDOFNumbering() );
    } else {
        buth = std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
        buth->setValues( 0.0 );
        buth->build();
    }

    if ( _phys_problem->getMesh()->isParallel() )
        CALLO_AP_ASSEMBLY_VECTOR( buth->getName() );

    return buth;
};