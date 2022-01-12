/**
 * @file BaseMesh.cxx
 * @brief Implementation de BaseMesh
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "astercxx.h"

#include "aster_fort_mesh.h"
#include "aster_fort_superv.h"
#include "aster_fort_utils.h"

#include "Meshes/BaseMesh.h"
#include "PythonBindings/LogicalUnitManager.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/Exceptions.h"
#include "Supervis/ResultNaming.h"
#include "Utilities/CapyConvertibleValue.h"
#include "Utilities/Tools.h"


int BaseMesh::getNumberOfNodes() const {
    if ( isEmpty() )
        return 0;
    if ( !_dimensionInformations->updateValuePointer() )
        return 0;
    return ( *_dimensionInformations )[0];
}

int BaseMesh::getNumberOfCells() const {
    if ( isEmpty() )
        return 0;
    if ( !_dimensionInformations->updateValuePointer() )
        return 0;
    return ( *_dimensionInformations )[2];
}

int BaseMesh::getDimension() const {
    if ( isEmpty() )
        return 0;
    if ( !_dimensionInformations->updateValuePointer() )
        return 0;

    const int dimGeom = ( *_dimensionInformations )[5];

    if ( dimGeom == 3 ) {
        const std::string typeco( "MAILLAGE" );
        ASTERINTEGER repi = 0, ier = 0;
        JeveuxChar32 repk( " " );
        const std::string arret( "F" );
        const std::string questi( "DIM_GEOM" );

        CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );

        return repi;
    }
    return dimGeom;
}

bool BaseMesh::readMeshFile( const std::string &fileName, const std::string &format ) {
    FileType type = Ascii;
    if ( format == "MED" )
        type = Binary;
    LogicalUnitFile file1( fileName, type, Old );

    SyntaxMapContainer syntax;

    if ( format == "GIBI" || format == "GMSH" ) {
        // Fichier temporaire
        LogicalUnitFile file2( "", Ascii, New );
        std::string preCmd = "PRE_" + format;
        ASTERINTEGER op2 = 47;
        if ( format == "GIBI" )
            op2 = 49;

        CommandSyntax *cmdSt2 = new CommandSyntax( preCmd );
        SyntaxMapContainer syntax2;
        syntax2.container["UNITE_" + format] = file1.getLogicalUnit();
        syntax2.container["UNITE_MAILLAGE"] = file2.getLogicalUnit();
        cmdSt2->define( syntax2 );

        try {
            CALL_EXECOP( &op2 );
        } catch ( ... ) {
            throw;
        }
        delete cmdSt2;
        syntax.container["FORMAT"] = "ASTER";
        syntax.container["UNITE"] = file2.getLogicalUnit();

        CommandSyntax cmdSt( "LIRE_MAILLAGE" );
        cmdSt.setResult( ResultNaming::getCurrentName(), "MAILLAGE" );

        cmdSt.define( syntax );

        try {
            ASTERINTEGER op = 1;
            CALL_EXECOP( &op );
        } catch ( ... ) {
            throw;
        }
    } else {
        syntax.container["FORMAT"] = format;
        syntax.container["UNITE"] = file1.getLogicalUnit();

        CommandSyntax cmdSt( "LIRE_MAILLAGE" );
        cmdSt.setResult( ResultNaming::getCurrentName(), "MAILLAGE" );

        cmdSt.define( syntax );

        ASTERINTEGER op = 1;
        CALL_EXECOP( &op );
    }

    return build();
}

bool BaseMesh::readMedFile( const std::string &fileName ) {
    return readMeshFile( fileName, "MED" );
}

const JeveuxCollectionLong BaseMesh::getMedConnectivity() const {
    JeveuxChar24 objv( " " );
    CALLO_GET_MED_CONNECTIVITY( getName(), objv );
    JeveuxCollectionLong result( objv.toString() );
    return result;
}

const JeveuxVectorLong BaseMesh::getMedCellsTypes() const {
    JeveuxChar24 objv( " " );
    CALLO_GET_MED_TYPES( getName(), objv );
    JeveuxVectorLong result( objv.toString() );
    return result;
}

bool BaseMesh::printMedFile( const std::string fileName, bool local ) const
{
    LogicalUnitFile a( fileName, Binary, New );
    ASTERINTEGER retour = a.getLogicalUnit();
    CommandSyntax cmdSt( "IMPR_RESU" );

    SyntaxMapContainer dict;

    if( isParallel() || isConnection()  )
        dict.container["PROC0"] = "NON";
    else
        dict.container["PROC0"] = "OUI";

    dict.container["FORMAT"] = "MED";
    dict.container["UNITE"] = retour;

    ListSyntaxMapContainer listeResu;
    SyntaxMapContainer dict2;
    dict2.container["MAILLAGE"] = getName();
    listeResu.push_back( dict2 );
    dict.container["RESU"] = listeResu;

    if ( !local && isParallel() )
        dict.container["FICHIER_UNIQUE"] = "OUI";

    cmdSt.define( dict );

    try {
        ASTERINTEGER op = 39;
        CALL_EXECOP( &op );
    } catch ( ... ) {
        throw;
    }

    return true;
};

const JeveuxCollectionLong BaseMesh::getInverseConnectivity() const {
    JeveuxChar24 objv( DataStructureNaming::getNewName( 24 ) );
    std::string base( "G" );

    ASTERINTEGER listCell;
    ASTERINTEGER nbCell = 0;

    CALLO_CNCINV( getName(), &listCell, &nbCell, base, objv );
    JeveuxCollectionLong result( objv.toString() );
    return result;
};

std::string BaseMesh::getNodeName( const ASTERINTEGER& index) const
{
    return trim(_nameOfNodes->getStringFromIndex(index));
};

std::string BaseMesh::getCellName( const ASTERINTEGER& index) const
{
    return trim(_nameOfCells->getStringFromIndex(index));
};
