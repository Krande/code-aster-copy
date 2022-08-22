/**
 * @file DiscreteComputation.cxx
 * @brief Implementation of class DiscreteComputation
 * @section LICENCE
 *   Copyright (C) 1991 2022  EDF R&D                www.code-aster.org
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

#include "aster_fort_calcul.h"
#include "aster_fort_superv.h"
#include "astercxx.h"

#include "Discretization/Calcul.h"
#include "Discretization/DiscreteComputation.h"
#include "Loads/DirichletBC.h"
#include "Loads/MechanicalLoad.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Modeling/XfemModel.h"
#include "Utilities/Tools.h"

/** @brief Compute AFFE_CHAR_THER TEMP_IMPO */
bool DiscreteComputation::addTherImposedTerms( ElementaryVectorRealPtr elemVect,
                                               const ASTERDOUBLE time, const ASTERDOUBLE time_step,
                                               const ASTERDOUBLE theta ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    // Init
    ASTERINTEGER iload = 1;
    bool has_load = false;
    std::string load_option = "";

    // Setup
    const std::string calcul_option( "CHAR_THER" );
    elemVect->prepareCompute( calcul_option );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto listOfLoads = _phys_problem->getListOfLoads();

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );

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
                calcul->addTimeField( "PTEMPSR", time, time_step, theta );
                calcul->addOutputElementaryTerm( "PVECTTR",
                                                 std::make_shared< ElementaryTermReal >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                    elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                                 iload );
                    has_load = true;
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

    return has_load;
}

/** @brief Compute CHAR_THER_EVOL */
FieldOnNodesRealPtr
DiscreteComputation::transientThermalLoad( const ASTERDOUBLE time, const ASTERDOUBLE time_step,
                                           const ASTERDOUBLE theta,
                                           const FieldOnCellsRealPtr _externVarField,
                                           const FieldOnNodesRealPtr _previousNodalField ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );
    AS_ASSERT( _previousNodalField->exists() );

    auto elemVect = std::make_shared< ElementaryVectorReal >( _phys_problem );

    // Setup
    const std::string calcul_option( "CHAR_THER_EVOL" );
    elemVect->prepareCompute( calcul_option );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto listOfLoads = _phys_problem->getListOfLoads();

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );
    calcul->clearInputs();
    calcul->clearOutputs();

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    calcul->addTimeField( "PTEMPSR", time, time_step, theta );
    calcul->addInputField( "PTEMPER", _previousNodalField );

    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
    }
    if ( _externVarField ) {
        calcul->addInputField( "PVARCPR", _externVarField );
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

bool DiscreteComputation::addTherNeumannTerms(
    ElementaryVectorRealPtr elemVect, const ASTERDOUBLE time, const ASTERDOUBLE time_step,
    const ASTERDOUBLE theta, const FieldOnCellsRealPtr _externVarField,
    const FieldOnNodesRealPtr _previousNodalField ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );
    AS_ASSERT( _previousNodalField->exists() );

    // Init
    ASTERINTEGER iload = 1;
    bool has_load = false;
    std::string load_option = "";

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
    calcul->clearInputs();
    calcul->clearOutputs();

    auto therLoadReal = listOfLoads->getThermalLoadsReal();
    for ( const auto &load : therLoadReal ) {
        auto load_FEDesc = load->getFiniteElementDescriptor();
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
            CALLO_RSINCH( evol_char_name, para_flun, access_var, &time,
                          evol_flow_xyz_field->getName(), extr_right, extr_left, &stop, base,
                          &iret );

            if ( iret >= 2 ) {

                FieldOnCellsRealPtr evol_exchange_field =
                    std::make_shared< FieldOnCellsReal >( model_FEDesc );
                FieldOnCellsRealPtr evol_ext_temp_field =
                    std::make_shared< FieldOnCellsReal >( model_FEDesc );

                CALLO_RSINCH( evol_char_name, para_coefh, access_var, &time,
                              evol_exchange_field->getName(), extr_right, extr_left, &stop, base,
                              &iret );

                if ( iret >= 2 ) {
                    AS_ABORT( "Cannot find COEF_H in EVOL_CHAR " + evol_char_name + " at time " +
                              std::to_string( time ) );
                }

                CALLO_RSINCH( evol_char_name, para_text, access_var, &time,
                              evol_ext_temp_field->getName(), extr_right, extr_left, &stop, base,
                              &iret );

                if ( iret >= 2 ) {
                    AS_ABORT( "Cannot find T_EXT in EVOL_CHAR " + evol_char_name + " at time " +
                              std::to_string( time ) );
                }

                calcul->setOption( "CHAR_THER_TEXT_R" );
                calcul->setFiniteElementDescriptor( model_FEDesc );
                calcul->clearInputs();
                calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
                calcul->addInputField( "PTEMPER", _previousNodalField );
                calcul->addTimeField( "PTEMPSR", time, time_step, theta );
                calcul->addInputField( "PCOEFHR", evol_exchange_field );
                calcul->addInputField( "PT_EXTR", evol_ext_temp_field );
                calcul->clearOutputs();
                calcul->addOutputElementaryTerm( "PVECTTR",
                                                 std::make_shared< ElementaryTermReal >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                    elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                                 iload );
                    has_load = true;
                }
            } else {
                calcul->setOption( "CHAR_THER_FLUN_R" );
                calcul->setFiniteElementDescriptor( model_FEDesc );
                calcul->clearInputs();
                calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
                calcul->addInputField( "PFLUXNR", evol_flow_xyz_field );
                calcul->clearOutputs();
                calcul->addOutputElementaryTerm( "PVECTTR",
                                                 std::make_shared< ElementaryTermReal >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                    elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                                 iload );
                    has_load = true;
                }
            }
        }

        // Termes ECHANGE
        if ( load->hasLoadField( "COEFH" ) && load->hasLoadField( "T_EXT" ) ) {
            auto exchange_field = load->getConstantLoadField( "COEFH" );
            auto ext_temp_field = load->getConstantLoadField( "T_EXT" );

            calcul->setOption( "CHAR_THER_TEXT_R" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            if ( _previousNodalField ) {
                calcul->addInputField( "PTEMPER", _previousNodalField );
            }
            calcul->addTimeField( "PTEMPSR", time, time_step, theta );
            calcul->addInputField( "PCOEFHR", exchange_field );
            calcul->addInputField( "PT_EXTR", ext_temp_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
                has_load = true;
            }
        }

        // Termes FLUX XYZ
        if ( load->hasLoadField( "FLURE" ) ) {
            auto flow_xyz_field = load->getConstantLoadField( "FLURE" );
            calcul->setOption( "CHAR_THER_FLUN_R" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addInputField( "PFLUXNR", flow_xyz_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
                has_load = true;
            }
        }

        // Termes FLUX NORM
        if ( load->hasLoadField( "FLUR2" ) ) {
            auto flow_nor_field = load->getConstantLoadField( "FLUR2" );
            calcul->setOption( "CHAR_THER_FLUX_R" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addTimeField( "PTEMPSR", time, time_step, theta );
            calcul->addInputField( "PFLUXVR", flow_nor_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
                has_load = true;
            }
        }

        // Termes SOURCE
        if ( load->hasLoadField( "SOURE" ) ) {
            auto source_field = load->getConstantLoadField( "SOURE" );
            calcul->setOption( "CHAR_THER_SOUR_R" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addTimeField( "PTEMPSR", time, time_step, theta );
            if ( _externVarField ) {
                calcul->addInputField( "PVARCPR", _externVarField );
            }
            calcul->addInputField( "PSOURCR", source_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
                has_load = true;
            }
        }
        // Termes SOURCE CALCULEE
        if ( load->hasLoadField( "SOURC" ) ) {
            auto computed_source_field = load->getLoadField( "SOURC" );
            calcul->setOption( "CHAR_THER_SOUR_R" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addTimeField( "PTEMPSR", time, time_step, theta );
            if ( _externVarField ) {
                calcul->addInputField( "PVARCPR", _externVarField );
            }
            calcul->addInputField( "PSOURCR", computed_source_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
                has_load = true;
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
            if ( _previousNodalField ) {
                calcul->addInputField( "PTEMPER", _previousNodalField );
            }
            calcul->addTimeField( "PTEMPSR", time, time_step, theta );
            calcul->addInputField( "PHECHPR", wall_exchange_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
                has_load = true;
            }
        }

        // Termes PRE_GRAD_TEMP
        if ( load->hasLoadField( "GRAIN" ) ) {
            auto pregrad_field = load->getConstantLoadField( "GRAIN" );
            calcul->setOption( "CHAR_THER_GRAI_R" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            if ( _externVarField ) {
                calcul->addInputField( "PVARCPR", _externVarField );
            }
            if ( currMater ) {
                calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
            }
            calcul->addInputField( "PGRAINR", pregrad_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
                has_load = true;
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
            if ( _previousNodalField ) {
                calcul->addInputField( "PTEMPER", _previousNodalField );
            }
            calcul->addTimeField( "PTEMPSR", time, time_step, theta );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->addInputField( "PCOEFHF", exchange_field );
            calcul->addInputField( "PT_EXTF", ext_temp_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
                has_load = true;
            }
        }

        // Termes FLUX XYZ
        if ( load->hasLoadField( "FLURE" ) ) {
            auto flow_xyz_field = load->getConstantLoadField( "FLURE" );
            calcul->setOption( "CHAR_THER_FLUN_F" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addTimeField( "PTEMPSR", time, time_step, theta );
            calcul->addInputField( "PFLUXNF", flow_xyz_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
                has_load = true;
            }
        }

        // Termes FLUX NORM
        if ( load->hasLoadField( "FLUR2" ) ) {
            auto flow_nor_field = load->getConstantLoadField( "FLUR2" );
            calcul->setOption( "CHAR_THER_FLUX_F" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addTimeField( "PTEMPSR", time, time_step, theta );
            calcul->addInputField( "PFLUXVF", flow_nor_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
                has_load = true;
            }
        }

        if ( load->hasLoadField( "SOURE" ) ) {
            auto source_field = load->getConstantLoadField( "SOURE" );
            calcul->setOption( "CHAR_THER_SOUR_F" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addTimeField( "PTEMPSR", time, time_step, theta );
            if ( _externVarField ) {
                calcul->addInputField( "PVARCPR", _externVarField );
            }
            calcul->addInputField( "PSOURCF", source_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
                has_load = true;
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
            if ( _previousNodalField ) {
                calcul->addInputField( "PTEMPER", _previousNodalField );
            }
            calcul->addTimeField( "PTEMPSR", time, time_step, theta );
            calcul->addInputField( "PHECHPF", wall_exchange_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
                has_load = true;
            }
        }

        // Termes PRE_GRAD_TEMP
        if ( load->hasLoadField( "GRAIN" ) ) {
            auto pregrad_field = load->getConstantLoadField( "GRAIN" );
            calcul->setOption( "CHAR_THER_GRAI_F" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addTimeField( "PTEMPSR", time, time_step, theta );
            if ( _externVarField ) {
                calcul->addInputField( "PVARCPR", _externVarField );
            }
            if ( currMater ) {
                calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
            }
            calcul->addInputField( "PGRAINF", pregrad_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
                has_load = true;
            }
        }

        iload++;
    }

    return has_load;
}
