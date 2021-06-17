/**
 * @file MaterialFieldBuilder.cxx
 * @brief Implementation de MaterialFieldBuilder::build
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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
#include <stdexcept>
#include <typeinfo>

#include "Materials/MaterialFieldBuilder.h"
#include "Supervis/CommandSyntax.h"
#include "Utilities/SyntaxDictionary.h"

void MaterialFieldBuilder::buildClass( MaterialField &curMater,
                                    const ListOfExternalStateVariablesPtr &curExternalStateVariable)
{
    SyntaxMapContainer dict;

    if ( !curMater._mesh )
        throw std::runtime_error( "Mesh is undefined" );
    dict.container["MAILLAGE"] = curMater._mesh->getName();
    if( curMater._model )
        dict.container["MODELE"] = curMater._model->getName();

    ListSyntaxMapContainer listeAFFE;
    for ( auto &curIter : curMater._materialsFieldEntity ) {
        SyntaxMapContainer dict2;
        VectorString listOfMater;
        for ( const auto &curIter2 : curIter.first )
            listOfMater.push_back( curIter2->getName() );
        dict2.container["MATER"] = listOfMater;
        const MeshEntityPtr &tmp = curIter.second;
        if ( tmp->getType() == AllMeshEntitiesType )
            dict2.container["TOUT"] = "OUI";
        else if ( tmp->getType() == GroupOfCellsType )
            dict2.container["GROUP_MA"] = ( curIter.second )->getNames();
        else if ( tmp->getType() == GroupOfNodesType )
            dict2.container["GROUP_NO"] = ( curIter.second )->getNames();
        else if ( tmp->getType() == CellType )
            dict2.container["MAILLE"] = ( curIter.second )->getNames();
        else
            throw std::runtime_error( "Support entity undefined" );
        listeAFFE.push_back( dict2 );
    }
    dict.container["AFFE"] = listeAFFE;

    ListSyntaxMapContainer listeAFFE_COMPOR;
    for ( auto &curIter : curMater._behaviours ) {
        SyntaxMapContainer dict2;
        dict2.container["COMPOR"] = curIter.first->getName();
        const MeshEntityPtr &tmp = curIter.second;
        if ( tmp->getType() == AllMeshEntitiesType )
            dict2.container["TOUT"] = "OUI";
        else if ( tmp->getType() == GroupOfCellsType )
            dict2.container["GROUP_MA"] = ( curIter.second )->getName();
        else if ( tmp->getType() == CellType )
            dict2.container["MAILLE"] = ( curIter.second )->getName();
        else
            throw std::runtime_error( "Support entity undefined" );
        listeAFFE_COMPOR.push_back( dict2 );
    }
    dict.container["AFFE_COMPOR"] = listeAFFE_COMPOR;

    if( curExternalStateVariable != nullptr ) {
        ListSyntaxMapContainer listeAFFE_VARC;
        for ( auto &curIter : curExternalStateVariable->_externalVars ) {
            SyntaxMapContainer dict2;

            const auto &externalVar = ( *curIter.first );
            dict2.container["NOM_VARC"] = externalVar.getVariableName();
            const auto &inputField = externalVar.getValue();
            const auto &evolParam = externalVar.getEvolutionParameter();
            if ( inputField != nullptr )
                dict2.container["CHAM_GD"] = inputField->getName();

            if ( evolParam != nullptr ) {
                dict2.container["EVOL"] = evolParam->getTransientResult()->getName();
                dict2.container["PROL_DROITE"] = evolParam->getRightExtension();
                dict2.container["PROL_GAUCHE"] = evolParam->getLeftExtension();
                if ( evolParam->getFieldName() != "" )
                    dict2.container["NOM_CHAM"] = evolParam->getFieldName();
                if ( evolParam->getTimeFormula() != nullptr )
                    dict2.container["FONC_INST"] = evolParam->getTimeFormula()->getName();
                if ( evolParam->getTimeFunction() != nullptr )
                    dict2.container["FONC_INST"] = evolParam->getTimeFunction()->getName();
            }
            if ( externalVar.hasReferenceValue() )
                dict2.container["VALE_REF"] = externalVar.getReferenceValue();

            const MeshEntityPtr &tmp = curIter.second;
            if ( tmp->getType() == AllMeshEntitiesType )
                dict2.container["TOUT"] = "OUI";
            else if ( tmp->getType() == GroupOfCellsType )
                dict2.container["GROUP_MA"] = ( curIter.second )->getName();
            else if ( tmp->getType() == CellType )
                dict2.container["MAILLE"] = ( curIter.second )->getName();
            else
                throw std::runtime_error( "Support entity undefined" );
            listeAFFE_VARC.push_back( dict2 );
        }
        dict.container["AFFE_VARC"] = listeAFFE_VARC;
    }

    auto syntax = CommandSyntax( "AFFE_MATERIAU" );
    syntax.setResult( curMater.getName(), curMater.getType() );
    syntax.define( dict );
    // Maintenant que le fichier de commande est pret, on appelle OP006
    try {
        ASTERINTEGER op = 6;
        CALL_EXECOP( &op );
    } catch ( ... ) {
        throw;
    }
};

MaterialFieldPtr
MaterialFieldBuilder::build( MaterialFieldPtr &curMater,
                            const ListOfExternalStateVariablesPtr &curExternalStateVariable)

{
    MaterialFieldBuilder::buildClass( *curMater, curExternalStateVariable);
    return curMater;
};
