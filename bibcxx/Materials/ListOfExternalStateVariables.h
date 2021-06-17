#ifndef EXTERNALVARIABLEFIELD_H_
#define EXTERNALVARIABLEFIELD_H_

/**
 * @file ExternalStateVariableField.h
 * @brief Fichier entete de la classe ExternalStateVariable
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

#include "Materials/BaseExternalStateVariables.h"
#include "astercxx.h"


/**
 * @class ListOfExternalStateVariables
 * @brief Input variable on mesh
 * @author Nicolas Sellenet
 */
class ListOfExternalStateVariables {
    friend class MaterialFieldBuilder;

  private:
    /** @typedef std::list d'une std::pair de MeshEntityPtr */
    typedef std::vector< std::pair< BaseExternalStateVariablePtr, MeshEntityPtr > >
        VectorOfexternalVarAndGrps;
    /** @typedef Definition de la valeur contenue dans un VectorOfexternalVarAndGrps */
    typedef VectorOfexternalVarAndGrps::value_type VectorOfexternalVarAndGrpsValue;
    /** @typedef Definition d'un iterateur sur VectorOfexternalVarAndGrps */
    typedef VectorOfexternalVarAndGrps::iterator VectorOfexternalVarAndGrpsIter;

    /** @brief Vector of BaseExternalStateVariable */
    VectorOfexternalVarAndGrps _externalVars;
    /** @brief Maillage sur lequel repose la sd_cham_mater */
    BaseMeshPtr _mesh;

  public:
    ListOfExternalStateVariables( const MeshPtr &mesh ) : _mesh( mesh ){};

    ListOfExternalStateVariables( const SkeletonPtr &mesh ) : _mesh( mesh ){};

#ifdef ASTER_HAVE_MPI
    ListOfExternalStateVariables( const ParallelMeshPtr &mesh ) : _mesh( mesh ){};
#endif /* ASTER_HAVE_MPI */

    /**
     * @brief Add an input variable on all mesh
     */
    template < class ExternalStateVariablePtr >
    void addExternalStateVariableOnMesh( const ExternalStateVariablePtr &curBehav ) {
        _externalVars.push_back(
            VectorOfexternalVarAndGrpsValue( curBehav, MeshEntityPtr( new AllMeshEntities() ) ) );
    };

    /**
     * @brief Add an input variable on a group of mesh
     */
    template < class ExternalStateVariablePtr >
    void addExternalStateVariableOnGroupOfCells( const ExternalStateVariablePtr &curBehav,
                                            const std::string &nameOfGroup ) {
        if ( !_mesh )
            throw std::runtime_error( "Mesh is not defined" );
        if ( !_mesh->hasGroupOfCells( nameOfGroup ) )
            throw std::runtime_error( nameOfGroup + "not in mesh" );

        _externalVars.push_back( VectorOfexternalVarAndGrpsValue(
            curBehav, MeshEntityPtr( new GroupOfCells( nameOfGroup ) ) ) );
    };

    /**
     * @brief Add an input variable on an element
     */
    template < class ExternalStateVariablePtr >
    void addExternalStateVariableOnCell( const ExternalStateVariablePtr &curBehav,
                                    const std::string &nameOfCell ) {
        if ( !_mesh )
            throw std::runtime_error( "Mesh is not defined" );

        _externalVars.push_back(
            VectorOfexternalVarAndGrpsValue( curBehav, MeshEntityPtr( new Cell( nameOfCell ) ) ) );
    };
};

/**
 * @typedef ListOfExternalStateVariablesPtr
 * @brief Pointeur intelligent vers un ListOfExternalStateVariables
 */
typedef boost::shared_ptr< ListOfExternalStateVariables > ListOfExternalStateVariablesPtr;

#endif /* EXTERNALVARIABLEFIELD_H_ */
