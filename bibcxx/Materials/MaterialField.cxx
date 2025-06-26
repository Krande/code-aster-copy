/**
 * @file MaterialField.cxx
 * @brief Implementation de MaterialField
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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

#include "Materials/MaterialField.h"

#include "aster_fort_superv.h"

#include "Results/TransientResult.h"
#include "Supervis/CommandSyntax.h"
#include "Utilities/SyntaxDictionary.h"

MaterialField::MaterialField( const std::string &name, const BaseMeshPtr &mesh )
    : _mesh( mesh ),
      _model( nullptr ),
      DataStructure( name, 8, "CHAM_MATER" ),
      _champ_mat( ConstantFieldOnCellsChar8Ptr(
          new ConstantFieldOnCellsChar8( getName() + ".CHAMP_MAT ", mesh ) ) ),
      _compor( ConstantFieldOnCellsRealPtr(
          new ConstantFieldOnCellsReal( getName() + ".COMPOR ", mesh ) ) ),
      _cvrcNom( JeveuxVectorChar8( getName() + ".CVRCNOM" ) ),
      _cvrcGd( JeveuxVectorChar8( getName() + ".CVRCGD" ) ),
      _cvrcVarc( JeveuxVectorChar8( getName() + ".CVRCVARC" ) ),
      _cvrcCmp( JeveuxVectorChar8( getName() + ".CVRCCMP" ) ) {};

MaterialField::MaterialField( const BaseMeshPtr &mesh, MaterialFieldPtr mater )
    : MaterialField( ResultNaming::getNewResultName(), mesh ) {
    _materialsOnMeshEntities = mater->_materialsOnMeshEntities;
    _behaviourOnMeshEntities = mater->_behaviourOnMeshEntities;
    _extStateVariablesOnMeshEntities = mater->_extStateVariablesOnMeshEntities;
};

listOfMaterials MaterialField::getVectorOfMaterial() const {
    listOfMaterials toReturn;
    for ( const auto &curIter : _materialsOnMeshEntities )
        for ( const auto &curIter2 : curIter.first )
            toReturn.push_back( curIter2 );
    return toReturn;
};

listOfPartOfMaterialField MaterialField::getVectorOfPartOfMaterialField() const {
    listOfPartOfMaterialField toReturn;
    for ( const auto &curIter : _materialsOnMeshEntities ) {
        PartOfMaterialFieldPtr toPush( new PartOfMaterialField( curIter.first, curIter.second ) );
        toReturn.push_back( toPush );
    }
    return toReturn;
};

ASTERINTEGER nameToId( const std::string &name ) {
    if ( name[0] == 'M' || name[1] == 'N' ) {
        const std::string tmp( name.substr( 1, name.size() - 1 ) );
        return std::atoi( tmp.c_str() );
    } else {
        return std::atoi( name.c_str() );
    }
};

MaterialPtr MaterialField::getMaterialOnCell( const std::string cellName ) const {
    // const auto cellNames = _mesh->getCellNameMap();
    ASTERINTEGER cellId = nameToId( cellName );
    _champ_mat->build();
    auto size = _champ_mat->size();
    ASTERINTEGER pos = size - 1;
    bool found = false;
    for ( ; pos >= 0; --pos ) {
        const auto zDesc = _champ_mat->getZoneDescription( pos );
        auto locType = zDesc.getLocalizationType();
        if ( locType == ConstantFieldOnZone::LocalizationType::AllMesh ) {
            found = true;
            break;
        } else if ( locType == ConstantFieldOnZone::LocalizationType::ListOfCells ) {
            auto cellsList = zDesc.getListOfCells();
            for ( const auto &id : cellsList ) {
                if ( id == cellId ) {
                    found = true;
                    break;
                }
            }
            if ( found )
                break;
        } else {
            throw std::runtime_error( "Localization type not allowed in MaterialField" );
        }
    }
    if ( !found )
        throw std::runtime_error( "Cell " + cellName + " not found in Material" );

    const auto valAndComp = _champ_mat->getValues( pos );
    const auto &val = valAndComp.getValues();
    const auto materName = val[0];
    MaterialPtr toReturn;
    found = false;
    auto tmp = materName.toString();
    for ( const auto &pairMatEnt : _materialsOnMeshEntities ) {
        const auto &curList = pairMatEnt.first;
        for ( const auto &mater : curList ) {
            if ( mater->getName() == tmp ) {
                toReturn = mater;
                found = true;
                break;
            }
        }
        if ( found )
            break;
    }
    return toReturn;
};

void MaterialField::addBehaviourOnMesh( BehaviourDefinitionPtr &curBehav ) {
    _behaviourOnMeshEntities.push_back(
        listOfBehavioursOnMeshValue( curBehav, MeshEntityPtr( new AllMeshEntities() ) ) );
}

void MaterialField::addBehaviourOnGroupOfCells( BehaviourDefinitionPtr &curBehav,
                                                VectorString namesOfGroup ) {
    if ( !_mesh )
        throw std::runtime_error( "Mesh is not defined" );
    for ( const auto &nameOfGroup : namesOfGroup )
        if ( !_mesh->hasGroupOfCells( nameOfGroup ) )
            throw std::runtime_error( nameOfGroup + " not in mesh" );

    _behaviourOnMeshEntities.push_back( listOfBehavioursOnMeshValue(
        curBehav, MeshEntityPtr( new GroupOfCells( namesOfGroup ) ) ) );
}

void MaterialField::addMultipleMaterialOnMesh( std::vector< MaterialPtr > curMaters ) {
    _materialsOnMeshEntities.push_back(
        listOfMaterialsOnMeshValue( curMaters, MeshEntityPtr( new AllMeshEntities() ) ) );
}

void MaterialField::addMaterialOnMesh( MaterialPtr &curMater ) {
    addMultipleMaterialOnMesh( std::vector< MaterialPtr >( { curMater } ) );
}

void MaterialField::addMultipleMaterialOnGroupOfCells( std::vector< MaterialPtr > curMaters,
                                                       VectorString namesOfGroup ) {
    if ( !_mesh )
        throw std::runtime_error( "Mesh is not defined" );
    for ( const auto &nameOfGroup : namesOfGroup )
        if ( !_mesh->hasGroupOfCells( nameOfGroup ) )
            throw std::runtime_error( nameOfGroup + " not in mesh" );

    _materialsOnMeshEntities.push_back( listOfMaterialsOnMeshValue(
        curMaters, MeshEntityPtr( new GroupOfCells( namesOfGroup ) ) ) );
}

void MaterialField::addMaterialOnGroupOfCells( MaterialPtr &curMater, VectorString namesOfGroup ) {
    addMultipleMaterialOnGroupOfCells( std::vector< MaterialPtr >( { curMater } ), namesOfGroup );
}

void MaterialField::addExternalStateVariable( ExternalStateVariablePtr &currExte ) {
    auto meshFromExte = currExte->getMesh();
    if ( meshFromExte != _mesh ) {
        throw std::runtime_error( "Mesh from external state variable is not the same" );
    }
    MeshEntityPtr meshEntity = currExte->getLocalization();
    if ( meshEntity->getType() == GroupOfCellsType ) {
        _extStateVariablesOnMeshEntities.push_back( listOfExternalVarOnMeshValue(
            currExte, MeshEntityPtr( new GroupOfCells( meshEntity->getName() ) ) ) );
    } else if ( meshEntity->getType() == AllMeshEntitiesType ) {
        _extStateVariablesOnMeshEntities.push_back(
            listOfExternalVarOnMeshValue( currExte, MeshEntityPtr( new AllMeshEntities() ) ) );
    } else {
        AS_ABORT( "Unknown entity" );
    }
}

bool MaterialField::hasExternalStateVariable( const std::string &name ) {
    if ( hasExternalStateVariable() ) {
        for ( auto curIter : getExtStateVariablesOnMeshEntities() ) {
            const auto externVar = curIter.first;
            if ( ExternalVariableTraits::getExternVarTypeStr( externVar->getType() ) == name ) {
                return true;
            }
        }
    } else {
        return false;
    }
    return false;
};

bool MaterialField::hasExternalStateVariable( const externVarEnumInt externVariType ) {
    if ( hasExternalStateVariable() ) {
        for ( auto curIter : getExtStateVariablesOnMeshEntities() ) {
            const auto externVar = curIter.first;
            if ( externVar->getType() == externVariType ) {
                return true;
            }
        }
    } else {
        return false;
    }
    return false;
};

bool MaterialField::hasExternalStateVariableForLoad() {
    if ( hasExternalStateVariable() ) {
        for ( auto curIter : getExtStateVariablesOnMeshEntities() ) {
            const auto externVar = curIter.first;
            if ( ExternalVariableTraits::externVarHasStrain( externVar->getType() ) ) {
                return true;
            }
        }
    } else {
        return false;
    }
    return false;
};

bool MaterialField::hasExternalStateVariableWithReference() {
    if ( hasExternalStateVariable() ) {
        for ( auto curIter : getExtStateVariablesOnMeshEntities() ) {
            const auto externVar = curIter.first;
            if ( ExternalVariableTraits::externVarHasRefeValue( externVar->getType() ) ) {
                return true;
            }
        }
    } else {
        return false;
    }
    return false;
};

ListSyntaxMapContainer MaterialField::syntaxForMaterial() {
    ListSyntaxMapContainer listAFFE;
    for ( auto &curIter : getMaterialsOnMeshEntities() ) {
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
        else
            throw std::runtime_error( "Support entity undefined" );
        listAFFE.push_back( dict2 );
    }
    return listAFFE;
}

ListSyntaxMapContainer MaterialField::syntaxForBehaviour() {
    ListSyntaxMapContainer listAFFE_COMPOR;
    for ( auto &curIter : getBehaviourOnMeshEntities() ) {
        SyntaxMapContainer dict2;
        dict2.container["COMPOR"] = curIter.first->getName();
        const MeshEntityPtr &tmp = curIter.second;
        if ( tmp->getType() == AllMeshEntitiesType )
            dict2.container["TOUT"] = "OUI";
        else if ( tmp->getType() == GroupOfCellsType )
            dict2.container["GROUP_MA"] = ( curIter.second )->getNames();
        else
            throw std::runtime_error( "Support entity undefined" );
        listAFFE_COMPOR.push_back( dict2 );
    }
    return listAFFE_COMPOR;
}

ListSyntaxMapContainer MaterialField::syntaxForExtStateVariables() {
    ListSyntaxMapContainer listAFFE_VARC;
    for ( auto &curIter : getExtStateVariablesOnMeshEntities() ) {
        SyntaxMapContainer dict2;
        const auto &externVar = ( *curIter.first );
        const auto externVariType = externVar.getType();
        dict2.container["NOM_VARC"] = ExternalVariableTraits::getExternVarTypeStr( externVariType );
        const auto &inputField = externVar.getField();
        const auto &evolParam = externVar.getEvolutionParameter();
        if ( inputField != nullptr ) {
            dict2.container["CHAM_GD"] = inputField->getName();
        };
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
        };
        if ( externVar.isSetRefe() ) {
            AS_ASSERT( ExternalVariableTraits::externVarHasRefeValue( externVariType ) );
            dict2.container["VALE_REF"] = externVar.getReferenceValue();
        } else {
            AS_ASSERT( !ExternalVariableTraits::externVarHasRefeValue( externVariType ) );
        }
        const MeshEntityPtr &meshEntity = curIter.second;
        if ( meshEntity->getType() == AllMeshEntitiesType ) {
            dict2.container["TOUT"] = "OUI";
        } else if ( meshEntity->getType() == GroupOfCellsType ) {
            dict2.container["GROUP_MA"] = ( curIter.second )->getName();
        } else {
            AS_ABORT( "Development error" );
        }
        listAFFE_VARC.push_back( dict2 );
    };
    return listAFFE_VARC;
}

bool MaterialField::build() {
    // Dictionnary for syntax
    SyntaxMapContainer commandSyntax;

    // Add model or mesh in dictionnary
    if ( !_mesh )
        throw std::runtime_error( "Mesh is undefined" );
    commandSyntax.container["MAILLAGE"] = _mesh->getName();
    if ( _model )
        commandSyntax.container["MODELE"] = _model->getName();

    // Add list of material parameters in dictionnary
    auto listAFFE = syntaxForMaterial();
    commandSyntax.container["AFFE"] = listAFFE;

    // Add list of multi-material parameters in dictionnary
    auto listAFFE_COMPOR = syntaxForBehaviour();
    commandSyntax.container["AFFE_COMPOR"] = listAFFE_COMPOR;

    // Add external state variables in dictionnary
    if ( hasExternalStateVariable() ) {
        auto listAFFE_VARC = syntaxForExtStateVariables();
        commandSyntax.container["AFFE_VARC"] = listAFFE_VARC;
    }

    // Final command
    auto commandFortran = CommandSyntax( "AFFE_MATERIAU" );
    commandFortran.setResult( getName(), getType() );
    commandFortran.define( commandSyntax );

    // Call Fortran command
    ASTERINTEGER op = 6;
    CALL_EXECOP( &op );

    // Generate C++ objects from Fortran objects for external state variables
    updateExtStateVariablesObjects();

    return true;
};

void MaterialField::updateExtStateVariablesObjects() {
    if ( _cvrcVarc->exists() ) {
        _cvrcVarc->updateValuePointer();
        ASTERINTEGER size = _cvrcVarc->size();
        for ( ASTERINTEGER i = 0; i < size; ++i ) {
            std::string varc_name = ( *_cvrcVarc )[i].toString();
            if ( _mapCvrcCard1.find( varc_name ) == _mapCvrcCard1.end() ) {
                _mapCvrcCard1[varc_name] = std::make_shared< ConstantFieldOnCellsReal >(
                    getName() + '.' + varc_name + ".1", _mesh );
                _mapCvrcCard2[varc_name] = std::make_shared< ConstantFieldOnCellsChar16 >(
                    getName() + '.' + varc_name + ".2", _mesh );
            }
        }
    }
};

bool MaterialField::updateInternalState() {
    updateExtStateVariablesObjects();
    return true;
};
