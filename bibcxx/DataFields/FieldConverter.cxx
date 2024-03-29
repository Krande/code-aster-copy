/**
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
#include "DataFields/FieldConverter.h"

FieldOnNodesReal toFieldOnNodes( const MeshCoordinatesField &field, const BaseMeshPtr mesh ) {
    FieldOnNodesReal chamno = FieldOnNodesReal();

    std::string type = "NOEU", base = "G";
    std::string prol = "NON", model = " ";

    CALLO_CHPCHD( field.getName(), type, mesh->getName(), prol, base, chamno.getName(), model );

    chamno.build( mesh );

    return chamno;
}

FieldOnNodesReal getRealPart( const FieldOnNodesComplex &field ) {

    auto newValues = JeveuxVectorReal( "&&TMP", field.size() );

    field.updateValuePointers();
    std::transform( field.getValues()->begin(), field.getValues()->end(), newValues.begin(),
                    []( ASTERCOMPLEX db ) { return db.real(); } );

    FieldOnNodesReal newField =
        FieldOnNodesReal( field.getEquationNumbering(), field.getReference(), newValues );

    return newField;
};

FieldOnNodesReal getImaginaryPart( const FieldOnNodesComplex &field ) {

    auto newValues = JeveuxVectorReal( "&&TMP", field.size() );

    field.updateValuePointers();
    std::transform( field.getValues()->begin(), field.getValues()->end(), newValues.begin(),
                    []( ASTERCOMPLEX db ) { return db.imag(); } );

    FieldOnNodesReal newField =
        FieldOnNodesReal( field.getEquationNumbering(), field.getReference(), newValues );

    return newField;
};

FieldOnNodesReal getRealPart( const FieldOnNodesReal &field ) { return field; };

FieldOnNodesReal getImaginaryPart( const FieldOnNodesReal &field ) {

    FieldOnNodesReal newField = FieldOnNodesReal( field );
    newField.setValues( 0. );

    return newField;
};
