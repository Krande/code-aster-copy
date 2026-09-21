/**
 * @file MedCouplingConversion.cxx
 * @brief Implementation de MedCouplingConversion
 * @author Francesco Bettonte
 * @section LICENCE
 *   Copyright (C) 1991 - 2026  EDF www.code-aster.org
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

#include "aster_pybind.h"

#include "Meshes/BaseMesh.h"

struct MedCouplingTypeInfo {
    int mc_type;
    int dim;
    int nb_nodes;
};

py::object getMedCouplingConversionData( const BaseMeshPtr &mesh ) {

    // Table de correspondance entre les types de mailles MED et leurs propriétés MEDCoupling.
    static const std::unordered_map< int, MedCouplingTypeInfo > med_to_mc = {
        { 1, { 0, 0, 1 } },     // POINT1
        { 102, { 1, 1, 2 } },   // SEG2
        { 103, { 2, 1, 3 } },   // SEG3
        { 104, { 10, 1, 4 } },  // SEG4
        { 203, { 3, 2, 3 } },   // TRI3
        { 206, { 6, 2, 6 } },   // TRI6
        { 207, { 7, 2, 7 } },   // TRI7
        { 204, { 4, 2, 4 } },   // QUAD4
        { 208, { 8, 2, 8 } },   // QUAD8
        { 209, { 9, 2, 9 } },   // QUAD9
        { 304, { 14, 3, 4 } },  // TETRA4
        { 310, { 20, 3, 10 } }, // TETRA10
        { 306, { 16, 3, 6 } },  // PENTA6
        { 315, { 25, 3, 15 } }, // PENTA15
        { 318, { 28, 3, 18 } }, // PENTA18
        { 305, { 15, 3, 5 } },  // PYRA5
        { 313, { 23, 3, 13 } }, // PYRA13
        { 308, { 18, 3, 8 } },  // HEXA8
        { 320, { 30, 3, 20 } }, // HEXA20
        { 327, { 27, 3, 27 } }  // HEXA27
    };

    // Structures de données pour stocker les informations du maillage triées par dimension.
    std::array< VectorLong, 4 > connectivity;
    std::array< VectorLong, 4 > connectivity_index;
    std::array< VectorLong, 4 > cells_renum;
    std::array< std::map< std::string, VectorLong >, 4 > groups_c;
    std::map< std::string, VectorLong > groups_n;

    JeveuxVectorLong cells_types = mesh->getMedCellsTypes();
    cells_types->updateValuePointer();
    const auto num_cells = cells_types->size();

    // Première passe sur les mailles pour collecter des informations globales.
    VectorInt cell_dim( num_cells );
    std::array< ASTERINTEGER, 4 > cells_per_dim_count = { 0, 0, 0, 0 };
    std::array< size_t, 4 > connectivity_total_size = { 0, 0, 0, 0 };

    for ( ASTERINTEGER i = 0; i < num_cells; ++i ) {
        const auto type_med = ( *cells_types )[i];
        if ( type_med == 0 ) {
            throw std::runtime_error( "Mesh contains non-med types and cannot be converted" );
        }
        const auto &type_info = med_to_mc.at( type_med );
        cell_dim[i] = type_info.dim;
        cells_per_dim_count[type_info.dim]++;
        connectivity_total_size[type_info.dim] += ( 1 + type_info.nb_nodes );
    }

    for ( int dim = 0; dim < 4; ++dim ) {
        connectivity[dim].reserve( connectivity_total_size[dim] );
        cells_renum[dim].reserve( cells_per_dim_count[dim] );
        connectivity_index[dim].reserve( cells_per_dim_count[dim] + 1 );
        // Premier element alloué ici
        connectivity_index[dim].push_back( 0 );
    }

    // Table pour mapper l'indice global d'une maille à son nouvel indice par dimension.
    VectorLong original_to_new_cell_idx( num_cells );
    std::array< ASTERINTEGER, 4 > current_cell_count_per_dim = { 0, 0, 0, 0 };

    JeveuxCollectionLong med_connectivity = mesh->getMedConnectivity();
    med_connectivity->build();

    // Seconde passe pour remplir les structures de connectivité.
    for ( ASTERINTEGER i = 0; i < num_cells; ++i ) {
        const auto type_med = ( *cells_types )[i];
        const auto &type_info = med_to_mc.at( type_med );
        const int dim = type_info.dim;

        // Ajoute le type de maille MEDCoupling.
        connectivity[dim].push_back( type_info.mc_type );

        // Ajoute les indices des nœuds (en base 0).
        auto nodes_med = ( *med_connectivity )[i + 1];
        nodes_med->updateValuePointer();
        for ( int j = 0; j < type_info.nb_nodes; ++j ) {
            connectivity[dim].push_back( ( *nodes_med )[j] - 1 );
        }

        // Calcule l'index de début de la maille suivante dans le tableau de connectivité.
        connectivity_index[dim].push_back( connectivity_index[dim].back() + 1 +
                                           type_info.nb_nodes );

        // Enregistre le nouvel indice de la maille dans la table de correspondance.
        original_to_new_cell_idx[i] = current_cell_count_per_dim[dim]++;

        cells_renum[dim].push_back( i + 1 );
    }

    // Traitement des groupes de mailles.
    for ( const auto &group_name : mesh->getGroupsOfCells() ) {
        const auto &cells_in_group = mesh->getCells( group_name );
        for ( const auto &cell_idx : cells_in_group ) {
            if ( cell_idx < 0 ) {
                throw std::runtime_error(
                    "Mesh contains non-med group id (negative) and cannot be converted" );
            }
            const int dim = cell_dim[cell_idx];
            const ASTERINTEGER new_idx = original_to_new_cell_idx[cell_idx];
            groups_c[dim][group_name].push_back( new_idx );
        }
    }

    // Traitement des groupes de nœuds.
    for ( const auto &group_name : mesh->getGroupsOfNodes() ) {
        const auto &nodes_in_group = mesh->getNodes( group_name, true );
        groups_n[group_name].reserve( nodes_in_group.size() );
        for ( const auto &node_idx : nodes_in_group ) {
            if ( node_idx < 0 ) {
                throw std::runtime_error(
                    "Mesh contains non-med group id (negative) and cannot be converted" );
            }
            groups_n[group_name].push_back( node_idx );
        }
    }

    // Conversion des données C++ en objets Python pour le retour.
    py::dict cells_dict;
    for ( int dim = 0; dim < 4; ++dim ) {
        if ( !connectivity[dim].empty() ) {
            // Crée un dictionnaire Python {dimension: (connectivité, index)}.
            cells_dict[py::int_( dim )] =
                py::make_tuple( py::cast( connectivity[dim] ), py::cast( connectivity_index[dim] ),
                                py::cast( cells_renum[dim] ) );
        }
    }

    py::dict groups_c_dict;
    for ( int dim = 0; dim < 4; ++dim ) {
        if ( !groups_c[dim].empty() ) {
            // Crée un dictionnaire Python {dimension: {nom_groupe: [indices]}}.
            groups_c_dict[py::int_( dim )] = py::cast( groups_c[dim] );
        }
    }

    // Convertit la map des groupes de nœuds en dictionnaire Python.
    py::dict groups_n_dict = py::cast( groups_n );

    // Retourne un tuple Python contenant les trois dictionnaires.
    return py::make_tuple( cells_dict, groups_c_dict, groups_n_dict );
}
