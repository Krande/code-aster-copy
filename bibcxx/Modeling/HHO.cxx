
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

#include "Modeling/HHO.h"

FieldOnNodesRealPtr HHO::projectOnLagrangeSpace( const FieldOnNodesRealPtr hho_field ) const {

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

FieldOnNodesRealPtr HHO::projectOnHHOSpace( const ASTERDOUBLE &value ) const {

    auto hho_field = std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
    hho_field->setValues( 0.0 );

    std::map< std::string, ASTERDOUBLE > map;

    auto model = _phys_problem->getModel();

    if ( model->isMechanical() ) {
        AS_ABORT( "TODO..." );
        // Comment faire pour les cells...
        map["HHO_U1"] = value;
        map["HHO_V1"] = value;
        map["HHO_W1"] = value;
    } else if ( model->isThermal() ) {
        map["HHO_C1"] = value;
        map["HHO_F1"] = value;
    } else {
        AS_ABORT( "Not implemented for HHO" );
    }

    hho_field->setValues( map );

    return hho_field;
};
