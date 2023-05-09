/**
 * @file SuperMesh.h
 * @brief Fichier entete de la classe SuperMesh
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

#include "Meshes/Mesh.h"

class DynamicMacroElement;
typedef std::shared_ptr< DynamicMacroElement > DynamicMacroElementPtr;
class StaticMacroElement;
typedef std::shared_ptr< StaticMacroElement > StaticMacroElementPtr;

class SuperMesh : public Mesh {

    /** @brief jeveux vector '.NOMACR' */
    JeveuxVectorLong _superElementName;
    /** @brief jeveux vector '.PARA_R' */
    JeveuxVectorReal _superElementPara;
    /** @brief jeveux vector '.SUPMAIL' */
    JeveuxCollectionLong _superElements;

    std::vector< DynamicMacroElementPtr > _dynamic_macro_elements;
    std::vector< StaticMacroElementPtr > _static_macro_elements;

  public:
    /**
     * @brief Constructor
     */
    SuperMesh();

    /**
     * @brief Constructor
     */
    SuperMesh( const std::string & );

    /**
     * @brief Add a DynamicMacroElement
     */
    bool addDynamicMacroElement( const DynamicMacroElementPtr & );

    /**
     * @brief Get all DynamicMacroElements
     */
    std::vector< DynamicMacroElementPtr > getDynamicMacroElements() const;

    /**
     * @brief Add a StaticMacroElement
     */
    bool addStaticMacroElement( const StaticMacroElementPtr & );

    /**
     * @brief Get all StaticMacroElements
     */
    std::vector< StaticMacroElementPtr > getStaticMacroElements() const;

    bool build();
};
