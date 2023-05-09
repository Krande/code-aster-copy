/**
 * @file SuperMesh.cxx
 * @brief Implementation de SuperMesh
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

#include "Meshes/SuperMesh.h"

#include "Modal/DynamicMacroElement.h"
#include "Modal/StaticMacroElement.h"
#include "Supervis/ResultNaming.h"

SuperMesh::SuperMesh() : SuperMesh( ResultNaming::getNewResultName() ) {};

SuperMesh::SuperMesh( const std::string &name )
    : Mesh( name ),
      _superElementName( JeveuxVectorLong( getName() + ".NOMACR" ) ),
      _superElementPara( JeveuxVectorReal( getName() + ".PARA_R" ) ),
      _superElements( JeveuxCollectionLong( getName() + ".SUPMAIL" ) ) {};

bool SuperMesh::addDynamicMacroElement( const DynamicMacroElementPtr &elem ) {
    _dynamic_macro_elements.push_back( elem );
    return true;
};

std::vector< DynamicMacroElementPtr > SuperMesh::getDynamicMacroElements() const {
    return _dynamic_macro_elements;
};

bool SuperMesh::addStaticMacroElement( const StaticMacroElementPtr &elem ) {
    _static_macro_elements.push_back( elem );
    return true;
};

std::vector< StaticMacroElementPtr > SuperMesh::getStaticMacroElements() const {
    return _static_macro_elements;
};

bool SuperMesh::build() {
    _superElements->build();
    return Mesh::build();
};
