/**
 * @file DirichletBC.cxx
 * @brief Implementation de DirichletBC
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

#include <typeinfo>
#include "astercxx.h"

#include "aster_fort_superv.h"
#include "Loads/DirichletBC.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/ResultNaming.h"

DirichletBCClass::DirichletBCClass( const std::string &type, const ModelPtr& model )
    : DataStructure( ResultNaming::getNewResultName(), 19, "CHAR_CINE" + type ),
      _intParam( JeveuxVectorLong( getName() + ".AFCI" ) ),
      _charParam( JeveuxVectorChar8( getName() + ".AFCK" ) ),
      _doubleParam( JeveuxVectorReal( getName() + ".AFCV" ) ), _isEmpty( true )
      {
          this->setModel(model);
      };

DirichletBCClass::DirichletBCClass(   const std::string &name,
                                            const std::string &type,
                                            const ModelPtr& model )
    : DataStructure( name, 19, "CHAR_CINE" + type ),
      _intParam( JeveuxVectorLong( getName() + ".AFCI" ) ),
      _charParam( JeveuxVectorChar8( getName() + ".AFCK" ) ),
      _doubleParam( JeveuxVectorReal( getName() + ".AFCV" ) ), _isEmpty( true )
      {
          this->setModel(model);
      };

bool DirichletBCClass::build() {
    std::string cmd = "AFFE_CHAR_CINE";
    if ( _listOfFunctionImposedTemperature.size() != 0 ) {
        cmd += "_F";
        throw std::runtime_error("Not implemented");
    }
    CommandSyntax cmdSt( cmd );
    cmdSt.setResult( ResultNaming::getCurrentName(), getType() );

    SyntaxMapContainer dict;
    if ( !_model )
        throw std::runtime_error( "Model is undefined" );
    dict.container["MODELE"] = _model->getName();

    // Definition de mot cle facteur MECA_IMPO
    if ( _listOfRealImposedDisplacement.size() != 0 ) {
        ListSyntaxMapContainer listeMecaImpo;
        for ( ListDispRealIter curIter = _listOfRealImposedDisplacement.begin();
              curIter != _listOfRealImposedDisplacement.end(); ++curIter ) {
            SyntaxMapContainer dict2;
            const MeshEntityPtr &tmp = curIter->getMeshEntityPtr();
            if ( tmp->getType() == AllMeshEntitiesType ) {
                dict2.container["TOUT"] = "OUI";
            } else {
                if ( tmp->getType() == GroupOfNodesType )
                    dict2.container["GROUP_NO"] = tmp->getName();
                else if ( tmp->getType() == GroupOfCellsType )
                    dict2.container["GROUP_MA"] = tmp->getName();
            }

            const std::string nomComp = curIter->getAsterCoordinateName();
            dict2.container[nomComp] = curIter->getValue();

            listeMecaImpo.push_back( dict2 );
        }

        dict.container["MECA_IMPO"] = listeMecaImpo;
    }
    // Definition de mot cle facteur THER_IMPO
    if ( _listOfRealImposedTemperature.size() != 0 ) {
        ListSyntaxMapContainer listeTempImpo;
        for ( ListDispTempIter curIter = _listOfRealImposedTemperature.begin();
              curIter != _listOfRealImposedTemperature.end(); ++curIter ) {
            SyntaxMapContainer dict2;
            const MeshEntityPtr &tmp = curIter->getMeshEntityPtr();
            if ( tmp->getType() == AllMeshEntitiesType ) {
                dict2.container["TOUT"] = "OUI";
            } else {
                if ( tmp->getType() == GroupOfNodesType )
                    dict2.container["GROUP_NO"] = tmp->getName();
                else if ( tmp->getType() == GroupOfCellsType )
                    dict2.container["GROUP_MA"] = tmp->getName();
            }
            const std::string nomComp = curIter->getAsterCoordinateName();
            dict2.container[nomComp] = curIter->getValue();

            listeTempImpo.push_back( dict2 );
        }

        dict.container["THER_IMPO"] = listeTempImpo;
    }
    cmdSt.define( dict );

    try {
        ASTERINTEGER op = 101;
        CALL_EXECOP( &op );
    } catch ( ... ) {
        throw;
    }
    _isEmpty = false;

    return true;
};
