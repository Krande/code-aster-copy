/**
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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

#include "DataFields/FieldConverter.h"
#include "Discretization/Calcul.h"

ModelPtr HHO::getModel() const {
    if ( _phys_problem ) {
        return _phys_problem->getModel();
    }

    return nullptr;
};

FunctionPtr HHO::_createFunc( const ASTERDOUBLE &value ) const {
    auto funct = std::make_shared< Function >();
    funct->setValues( { 1. }, { value } );
    funct->setResultName( "TOUTRESU" );
    funct->setParameterName( "TOUTPARA" );
    funct->setInterpolation( "LIN LIN" );
    funct->setExtrapolation( "CC" );
    funct->setAsConstant();

    return funct;
};

FieldOnNodesRealPtr HHO::projectOnLagrangeSpace( const FieldOnNodesRealPtr hho_field ) const {

    std::string option, para_name_in, para_name_out;

    auto model = this->getModel();

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

    calcul->addHHOField( model->getHHOModel() );

    // Add output terms
    auto exitField = std::make_shared< FieldOnCellsReal >( model );
    calcul->addOutputField( para_name_out, exitField );

    // Compute
    if ( model->existsFiniteElement() ) {
        calcul->compute();
    }

    return toFieldOnNodes( exitField );
};

FieldOnNodesRealPtr HHO::_projectOnHHOSpace( bool faces, const GenericFunctionPtr fct,
                                             ASTERDOUBLE time ) const {

    const std::string option = "HHO_PROJ_THER";
    auto model = this->getModel();
    auto mesh = model->getMesh();

    AS_ASSERT( model->isThermal() );

    auto calcul = std::make_unique< Calcul >( option );
    calcul->setModel( model );

    std::map< bool, std::string > dic_faces { { true, "ALL" }, { false, "CELL" } };

    auto funcField = std::make_shared< ConstantFieldOnCellsChar8 >( mesh );
    const std::string physicalName( "NEUT_K8" );
    funcField->allocate( physicalName );
    ConstantFieldOnZone a( mesh );
    ConstantFieldValues< JeveuxChar8 > b( { "Z1", "Z2" }, { fct->getName(), dic_faces[faces] } );
    funcField->setValueOnZone( a, b );

    // Input fields
    calcul->addInputField( "PGEOMER", mesh->getCoordinates() );
    calcul->addInputField( "PFUNC_R", funcField );
    calcul->addHHOField( model->getHHOModel() );
    calcul->addTimeField( "PINSTPR", time );

    // Output fields
    auto hho_elno = std::make_shared< FieldOnCellsReal >( model );
    calcul->addOutputField( "PTEMP_R", hho_elno );

    // Compute
    if ( model->existsFiniteElement() ) {
        calcul->compute();
    };

    return toFieldOnNodes( hho_elno );
};

FieldOnNodesRealPtr HHO::projectOnHHOSpace( const GenericFunctionPtr fct, ASTERDOUBLE time ) const {

    return _projectOnHHOSpace( true, fct, time );
};

FieldOnNodesRealPtr HHO::_projectOnHHOSpace( bool faces,
                                             const std::vector< GenericFunctionPtr > fct,
                                             ASTERDOUBLE time ) const {

    const std::string option = "HHO_PROJ_MECA";
    auto model = this->getModel();
    auto mesh = model->getMesh();
    auto dimMesh = mesh->getDimension();

    std::map< bool, std::string > dic_faces { { true, "ALL" }, { false, "CELL" } };

    AS_ASSERT( model->isMechanical() );
    AS_ASSERT( fct.size() == dimMesh );

    auto calcul = std::make_unique< Calcul >( option );
    calcul->setModel( model );

    auto funcField = std::make_shared< ConstantFieldOnCellsChar8 >( mesh );
    const std::string physicalName( "NEUT_K8" );
    funcField->allocate( physicalName );
    ConstantFieldOnZone a( mesh );

    if ( dimMesh == 2 ) {
        ConstantFieldValues< JeveuxChar8 > b(
            { "Z1", "Z2", "Z3" }, { fct[0]->getName(), fct[1]->getName(), dic_faces[faces] } );
        funcField->setValueOnZone( a, b );
    } else {
        ConstantFieldValues< JeveuxChar8 > b(
            { "Z1", "Z2", "Z3", "Z4" },
            { fct[0]->getName(), fct[1]->getName(), fct[2]->getName(), dic_faces[faces] } );
        funcField->setValueOnZone( a, b );
    }

    // Input fields
    calcul->addInputField( "PGEOMER", mesh->getCoordinates() );
    calcul->addInputField( "PFUNC_R", funcField );
    calcul->addHHOField( model->getHHOModel() );
    calcul->addTimeField( "PINSTPR", time );

    // Output fields
    auto hho_elno = std::make_shared< FieldOnCellsReal >( model );
    calcul->addOutputField( "PDEPL_R", hho_elno );

    // Compute
    if ( model->existsFiniteElement() ) {
        calcul->compute();
    };

    return toFieldOnNodes( hho_elno );
};

FieldOnNodesRealPtr HHO::projectOnHHOSpace( const std::vector< GenericFunctionPtr > fct,
                                            ASTERDOUBLE time ) const {
    return _projectOnHHOSpace( true, fct, time );
};

FieldOnNodesRealPtr HHO::projectOnHHOSpace( const ASTERDOUBLE &value ) const {
    return projectOnHHOSpace( _createFunc( value ) );
};

FieldOnNodesRealPtr HHO::projectOnHHOSpace( const VectorReal &values ) const {
    std::vector< GenericFunctionPtr > fct;
    for ( auto &val : values ) {
        fct.push_back( _createFunc( val ) );
    }

    return projectOnHHOSpace( fct );
};

FieldOnNodesRealPtr HHO::projectOnHHOCellSpace( const GenericFunctionPtr fct,
                                                ASTERDOUBLE time ) const {
    return _projectOnHHOSpace( false, fct, time );
};

FieldOnNodesRealPtr HHO::projectOnHHOCellSpace( const std::vector< GenericFunctionPtr > fct,
                                                ASTERDOUBLE time ) const {

    return _projectOnHHOSpace( false, fct, time );
};

FieldOnNodesRealPtr HHO::projectOnHHOCellSpace( const ASTERDOUBLE &value ) const {
    return projectOnHHOCellSpace( _createFunc( value ) );
};

FieldOnNodesRealPtr HHO::projectOnHHOCellSpace( const VectorReal &values ) const {
    std::vector< GenericFunctionPtr > fct;
    for ( auto &val : values ) {
        fct.push_back( _createFunc( val ) );
    }

    return projectOnHHOCellSpace( fct );
};

FieldOnNodesRealPtr HHO::projectOnHHOSpace( const FieldOnNodesRealPtr h1_field ) const {
    auto model = this->getModel();
    auto mesh = model->getMesh();

    const std::string option = model->isThermal() ? "HHO_PROJ2_THER" : "HHO_PROJ2_MECA";
    auto calcul = std::make_unique< Calcul >( option );
    calcul->setModel( model );

    // Input fields
    calcul->addInputField( "PGEOMER", mesh->getCoordinates() );
    calcul->addInputField( "PH1TP_R", h1_field );
    calcul->addHHOField( model->getHHOModel() );

    // Output fields
    auto hho_elno = std::make_shared< FieldOnCellsReal >( model );
    const std::string pname = model->isThermal() ? "PTEMP_R" : "PDEPL_R";
    calcul->addOutputField( pname, hho_elno );

    // Compute
    if ( model->existsFiniteElement() ) {
        calcul->compute();
    };

    return toFieldOnNodes( hho_elno );
};

FieldOnNodesRealPtr HHO::projectOnHHOCellSpace( const FieldOnCellsRealPtr field_elga ) const {
    auto model = this->getModel();
    auto mesh = model->getMesh();

    const std::string option = model->isThermal() ? "HHO_PROJ3_THER" : "HHO_PROJ3_MECA";
    auto calcul = std::make_unique< Calcul >( option );
    calcul->setModel( model );

    // Input fields
    calcul->addInputField( "PGEOMER", mesh->getCoordinates() );
    calcul->addInputField( "PQPTP_R", field_elga );
    calcul->addHHOField( model->getHHOModel() );

    // Output fields
    auto hho_elno = std::make_shared< FieldOnCellsReal >( model );
    const std::string pname = model->isThermal() ? "PTEMP_R" : "PDEPL_R";
    calcul->addOutputField( pname, hho_elno );

    // Compute
    if ( model->existsFiniteElement() ) {
        calcul->compute();
    };

    return toFieldOnNodes( hho_elno );
};

FieldOnCellsRealPtr HHO::evaluateAtQuadraturePoints( const FieldOnNodesRealPtr hho_field ) const {
    auto model = this->getModel();
    auto mesh = model->getMesh();

    const std::string option = model->isThermal() ? "TEMP_ELGA" : "DEPL_ELGA";
    auto calcul = std::make_unique< Calcul >( option );
    calcul->setModel( model );

    // Input fields
    calcul->addInputField( "PGEOMER", mesh->getCoordinates() );
    calcul->addInputField( "PTEMPER", hho_field );
    calcul->addHHOField( model->getHHOModel() );

    // Output fields
    auto hho_elga = std::make_shared< FieldOnCellsReal >( model );
    const std::string pname = model->isThermal() ? "PTEMP_R" : "PDEPL_R";
    calcul->addOutputField( pname, hho_elga );

    // Compute
    if ( model->existsFiniteElement() ) {
        calcul->compute();
    };

    hho_elga->updateValuePointers();

    return hho_elga;
};
