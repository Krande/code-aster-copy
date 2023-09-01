
/**
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

#include "PostProcessing/PostProcessing.h"

#include "DataFields/FieldOnCellsBuilder.h"
#include "Discretization/Calcul.h"
#include "Modeling/Model.h"

FieldOnCellsRealPtr PostProcessing::computeHydration( const FieldOnNodesRealPtr temp_prev,
                                                      const FieldOnNodesRealPtr temp_curr,
                                                      const ASTERDOUBLE time_prev,
                                                      const ASTERDOUBLE time_curr,
                                                      const FieldOnCellsRealPtr hydr_prev ) const {
    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    const std::string calcul_option( "HYDR_ELGA" );

    auto currModel = _phys_problem->getModel();
    auto currBehav = _phys_problem->getBehaviourProperty();
    auto currCodedMater = _phys_problem->getCodedMaterial();

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );

    // Add Input Field
    calcul->addInputField( "PCOMPOR", currBehav->getBehaviourField() );
    calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
    calcul->addTimeField( "PTEMPSR", time_curr, time_curr - time_prev, -1.0 );
    calcul->addInputField( "PTEMPMR", temp_prev );
    calcul->addInputField( "PTEMPPR", temp_curr );
    calcul->addInputField( "PHYDRMR", hydr_prev );

    // Create output fields
    auto hydr_curr = FieldOnCellsPtrBuilder< ASTERDOUBLE >( calcul->getFiniteElementDescriptor(),
                                                            "ELGA", "HYDR_R" );
    calcul->addOutputField( "PHYDRPR", hydr_curr );

    calcul->compute();

    return hydr_curr;
};
