/**
 *   Copyright (C) 1991 2025  EDF R&D                www.code-aster.org
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

#pragma once

#include "astercxx.h"

#include "aster_fort_ds.h"

#include "DataFields/ConstantFieldOnCells.h"
#include "DataFields/FieldOnCells.h"
#include "DataFields/FieldOnCellsBuilder.h"
#include "DataFields/FieldOnNodes.h"
#include "DataFields/MeshCoordinatesField.h"
#include "DataFields/SimpleFieldOnCells.h"
#include "DataFields/SimpleFieldOnNodes.h"
#include "Meshes/ConnectionMesh.h"
#include "Supervis/Exceptions.h"

//////// Convert to FieldOnNodes ////////////////
template < typename ValueType >
std::shared_ptr< FieldOnNodes< ValueType > >
toFieldOnNodes( const FieldOnCells< ValueType > &field ) {
    auto chamno = std::make_shared< FieldOnNodes< ValueType > >();

    std::string type = "NOEU", celmod = " ", base = "G";
    std::string prol = "OUI", model = " ";

    CALLO_CHPCHD( field.getName(), type, celmod, prol, base, chamno->getName(), model );

    chamno->build( field.getMesh() );

    return chamno;
};

template < typename ValueType >
std::shared_ptr< FieldOnNodes< ValueType > >
toFieldOnNodes( const std::shared_ptr< FieldOnCells< ValueType > > field ) {
    return toFieldOnNodes( *field );
};

template < typename ValueType >
std::shared_ptr< FieldOnNodes< ValueType > >
toFieldOnNodes( const SimpleFieldOnNodes< ValueType > &field ) {
    auto cham_no = std::make_shared< FieldOnNodes< ValueType > >();

    // Convert to CHAM_NO
    std::string prof = " ", prol0 = "NON", base = "G", kstop = "F";
    ASTERINTEGER iret = 0;
    CALLO_CNSCNO_WRAP( field.getName(), prof, prol0, base, cham_no->getName(), kstop, &iret );

    AS_ASSERT( iret == 0 );

    cham_no->build( field.getMesh() );
    return cham_no;
};

template < typename ValueType >
std::shared_ptr< FieldOnNodes< ValueType > >
toFieldOnNodes( const std::shared_ptr< SimpleFieldOnNodes< ValueType > > field ) {
    return toFieldOnNodes( *field );
}

template < typename ValueType >
std::shared_ptr< FieldOnNodes< ValueType > >
toFieldOnNodes( const SimpleFieldOnCells< ValueType > &field ) {
    return toFieldOnNodes( toSimpleFieldOnNodes( field ) );
};

template < typename ValueType >
std::shared_ptr< FieldOnNodes< ValueType > >
toFieldOnNodes( const std::shared_ptr< SimpleFieldOnCells< ValueType > > field ) {
    return toFieldOnNodes( *field );
};

////// Convert to SimpleFieldOnNodes ////////////
template < typename ValueType >
std::shared_ptr< SimpleFieldOnNodes< ValueType > >
toSimpleFieldOnNodes( const FieldOnNodes< ValueType > &field ) {
    auto toReturn = std::make_shared< SimpleFieldOnNodes< ValueType > >( field.getMesh() );
    const std::string resultName = toReturn->getName();
    const std::string inName = field.getName();
    CALLO_CNOCNS_WRAP( inName, JeveuxMemoryTypesNames[Permanent], resultName );
    toReturn->build();
    return toReturn;
};

template < typename ValueType >
std::shared_ptr< SimpleFieldOnNodes< ValueType > >
toSimpleFieldOnNodes( const std::shared_ptr< FieldOnNodes< ValueType > > field ) {
    return toSimpleFieldOnNodes( *field );
};

template < typename ValueType >
std::shared_ptr< SimpleFieldOnNodes< ValueType > >
toSimpleFieldOnNodes( const FieldOnCells< ValueType > &field ) {
    return toSimpleFieldOnNodes( toFieldOnNodes( field ) );
};

template < typename ValueType >
std::shared_ptr< SimpleFieldOnNodes< ValueType > >
toSimpleFieldOnNodes( const std::shared_ptr< FieldOnCells< ValueType > > field ) {
    return toSimpleFieldOnNodes( *field );
};

template < typename ValueType >
std::shared_ptr< SimpleFieldOnNodes< ValueType > >
toSimpleFieldOnNodes( const SimpleFieldOnCells< ValueType > &field ) {
    auto chs = std::make_shared< SimpleFieldOnNodes< ValueType > >( field.getMesh() );

    // Convert to CHAM_NO_S
    const std::string base = "G", kstop = "F";
    ASTERINTEGER iret = 0;
    std::string celpg = " ";

    if ( field.getLocalization() == "ELGA" ) {
        raiseAsterError( "Not implemented: conversion to ELGA" );
    }

    CALLO_CESCNS( field.getName(), celpg, base, chs->getName(), kstop, &iret );

    AS_ASSERT( iret == 0 );

    chs->build();
    return chs;
}

template < typename ValueType >
std::shared_ptr< SimpleFieldOnNodes< ValueType > >
toSimpleFieldOnNodes( const std::shared_ptr< SimpleFieldOnCells< ValueType > > field ) {
    return toSimpleFieldOnNodes( *field );
};

////// Convert to SimpleFieldOnCells ////////////
template < typename ValueType >
std::shared_ptr< SimpleFieldOnCells< ValueType > >
toSimpleFieldOnCells( const ConstantFieldOnCells< ValueType > &field,
                      const SimpleFieldOnCells< ValueType > &simpleFieldModel ) {
    auto chs = std::make_shared< SimpleFieldOnCells< ValueType > >( field.getMesh() );

    // Convert to CHAM_ELEM_S
    const std::string base = "G", kstop = "A", loc = "ELGA";
    ASTERINTEGER iret = 0;
    std::string cesmod = simpleFieldModel.getName();

    CALL_CARCES( field.getName(), loc, cesmod, base, chs->getName(), kstop, &iret );

    AS_ASSERT( iret == 0 );

    chs->build();
    return chs;
}

template < typename ValueType >
std::shared_ptr< SimpleFieldOnCells< ValueType > >
toSimpleFieldOnCells( const std::shared_ptr< ConstantFieldOnCells< ValueType > > field,
                      const SimpleFieldOnCells< ValueType > &simpleFieldModel ) {
    return toSimpleFieldOnCells( *field, simpleFieldModel );
}

template < typename ValueType >
std::shared_ptr< SimpleFieldOnCells< ValueType > >
toSimpleFieldOnCells( const FieldOnCells< ValueType > &field ) {
    auto toReturn = std::make_shared< SimpleFieldOnCells< ValueType > >( field.getMesh() );
    const std::string resultName = toReturn->getName();
    const std::string inName = field.getName();
    const std::string copyNan( "OUI" );
    CALLO_CELCES_WRAP( inName, JeveuxMemoryTypesNames[Permanent], resultName );
    toReturn->build();
    return toReturn;
}

template < typename ValueType >
std::shared_ptr< SimpleFieldOnCells< ValueType > >
toSimpleFieldOnCells( const std::shared_ptr< FieldOnCells< ValueType > > field ) {
    return toSimpleFieldOnCells( *field );
}

template < typename ValueType >
std::shared_ptr< FieldOnCells< ValueType > >
toFieldOnCells( const SimpleFieldOnCells< ValueType > &field, const FiniteElementDescriptorPtr fed,
                const std::string option = std::string(),
                const std::string nompar = std::string() ) {
    auto cham_elem = std::make_shared< FieldOnCells< ValueType > >();

    // Convert to CHAM_ELEM
    const std::string prol0 = "OUI", base = "G", kstop = "F";
    ASTERINTEGER iret = 0, nncp = 0;
    CALLO_CESCEL_WRAP( field.getName(), fed->getName(), option, nompar, prol0, &nncp, base,
                       cham_elem->getName(), kstop, &iret );

    AS_ASSERT( iret == 0 );

    cham_elem->build( { fed } );
    cham_elem->updateValuePointers();
    return cham_elem;
}

template < typename ValueType >
std::shared_ptr< FieldOnCells< ValueType > >
toFieldOnCells( const std::shared_ptr< FieldOnNodes< ValueType > > field,
                const FiniteElementDescriptorPtr fed, const std::string loc,
                const std::string option = std::string(),
                const std::string nompar = std::string() ) {
    return toFieldOnCells( *field, fed, loc, option, nompar );
}

template < typename ValueType >
std::shared_ptr< FieldOnCells< ValueType > >
toFieldOnCells( const FieldOnNodes< ValueType > &field, const FiniteElementDescriptorPtr fed,
                const std::string loc, const std::string option = std::string(),
                const std::string nompar = std::string() ) {
    auto cham_elem = std::make_shared< FieldOnCells< ValueType > >();

    // Convert to CHAM_ELEM
    std::string base = "G";
    std::string prol = "OUI", model = " ";

    auto data = FieldOnCellsPtrBuilder< ValueType >( fed, loc, field.getPhysicalQuantity() );

    CALLO_CHPCHD( field.getName(), loc, data->getName(), prol, base, cham_elem->getName(), model );

    cham_elem->build( { fed } );
    cham_elem->updateValuePointers();
    return cham_elem;
}

template < typename ValueType >
std::shared_ptr< FieldOnCells< ValueType > >
toFieldOnCells( const std::shared_ptr< SimpleFieldOnCells< ValueType > > field,
                const FiniteElementDescriptorPtr fed, const std::string option = std::string(),
                const std::string nompar = std::string() ) {
    return toFieldOnCells( *field, fed, option, nompar );
}

#ifdef ASTER_HAVE_MPI
FieldOnNodesRealPtr transferToConnectionMesh( const FieldOnNodesRealPtr, const ConnectionMeshPtr );
FieldOnNodesRealPtr transferFromConnectionToParallelMesh( const FieldOnNodesRealPtr,
                                                          const BaseMeshPtr );
#endif /* ASTER_HAVE_MPI */

using FieldOnNodesReal = FieldOnNodes< ASTERDOUBLE >;
using FieldOnNodesComplex = FieldOnNodes< ASTERCOMPLEX >;

FieldOnNodesReal toFieldOnNodes( const MeshCoordinatesField &field, const BaseMeshPtr mesh );

FieldOnNodesReal getRealPart( const FieldOnNodesComplex &field );
FieldOnNodesReal getImaginaryPart( const FieldOnNodesComplex &field );
FieldOnNodesReal getRealPart( const FieldOnNodesReal &field );
FieldOnNodesReal getImaginaryPart( const FieldOnNodesReal &field );
