
/**
 * @section LICENCE
 *   Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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

#include "PostProcessing/PostProcessing.h"

FieldOnNodesRealPtr PostProcessing::projectHHO( const FieldOnNodesRealPtr hho_field ) const {

    std::string option, para_name_in, para_name_out;

    auto model = _phys_problem->getModel();

    if ( model->isMechanical() ) {
        option = "HHO_DEPL_MECA";
        para_name_in = "PDEPLPR";
        para_name_out = "PDEPL_R";
    } else if ( model->isThermal() ) {
        option = "HHO_TEMP_THER";
        para_name_in = "PTMPCHF";
        para_name_out = "PTEMP_R";
    } else {
        AS_ABORT( "Not implemented for HHO" );
    }

    // Main object
    CalculPtr calcul = std::make_unique< Calcul >( option );
    calcul->setModel( model );

    // Add input fields
    calcul->addInputField( "PGEOMER", model->getMesh()->getCoordinates() );
    calcul->addInputField( para_name_in, hho_field );

    // Add output terms
    auto exitField = std::make_shared< FieldOnCellsReal >( model );
    calcul->addOutputField( para_name_out, exitField );

    // Compute
    if ( model->existsFiniteElement() ) {
        calcul->compute();
    }

    return exitField->toFieldOnNodes();
};
