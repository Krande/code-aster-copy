/**
 * @file MedCouplingConversion.cxx
 * @brief Implementation de MedCouplingConversion
 * @author Francesco Bettonte
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

/* person_in_charge: francesco.bettonte at edf.fr */

#include "aster_pybind.h"

#include "Meshes/BaseMesh.h"

py::object getMedCouplingConversionData( const BaseMeshPtr &mesh ) {

    std::map< int, VectorInt > med_to_mc;
    std::map< int, VectorLong > connectivity;
    std::map< int, VectorLong > connectivity_index;
    std::map< int, MapLong > corresponding_cells;
    std::map< int, std::map< std::string, VectorLong > > groups_c;
    std::map< std::string, VectorLong > groups_n;

    // TYPE_MED = { TYPE_MEDCOUPLING, DIM, NB_NODES }
    med_to_mc[1] = { 0, 0, 1 };     // POINT1
    med_to_mc[102] = { 1, 1, 2 };   // SEG2
    med_to_mc[103] = { 2, 1, 3 };   // SEG3
    med_to_mc[104] = { 10, 1, 4 };  // SEG4
    med_to_mc[203] = { 3, 2, 3 };   // TRI3
    med_to_mc[206] = { 6, 2, 6 };   // TRI6
    med_to_mc[207] = { 7, 2, 7 };   // TRI7
    med_to_mc[204] = { 4, 2, 4 };   // QUAD4
    med_to_mc[208] = { 8, 2, 8 };   // QUAD8
    med_to_mc[209] = { 9, 2, 9 };   // QUAD9
    med_to_mc[304] = { 14, 3, 4 };  // TETRA4
    med_to_mc[310] = { 20, 3, 10 }; // TETRA10
    med_to_mc[306] = { 16, 3, 6 };  // PENTA6
    med_to_mc[315] = { 25, 3, 15 }; // PENTA15
    med_to_mc[318] = { 28, 3, 18 }; // PENTA18
    med_to_mc[305] = { 15, 3, 5 };  // PYRA5
    med_to_mc[313] = { 23, 3, 13 }; // PYRA13
    med_to_mc[308] = { 18, 3, 8 };  // HEXA8
    med_to_mc[320] = { 30, 3, 20 }; // HEXA20
    med_to_mc[327] = { 27, 3, 27 }; // HEXA27

    JeveuxVectorLong cells_types = mesh->getMedCellsTypes();
    cells_types->updateValuePointer();

    for ( ASTERINTEGER i = 0; i < cells_types->size(); ++i ) {
        if ( 0 == ( *cells_types )[i] ) {
            throw std::runtime_error( "Mesh contains non-med types and cannot be converted" );
        }
    }

    JeveuxCollectionLong med_connectivity = mesh->getMedConnectivity();
    med_connectivity->build();

    // Tri des mailles par dimension
    for ( ASTERINTEGER i = 0; i < cells_types->size(); ++i ) {

        auto aster_index = i + 1;
        auto type_med = ( *cells_types )[i];
        auto nodes_med = ( *med_connectivity )[aster_index];
        nodes_med->updateValuePointer();

        int mc_type = med_to_mc[type_med][0];
        int dim = med_to_mc[type_med][1];
        int cell_size = med_to_mc[type_med][2];

        connectivity[dim].push_back( mc_type );
        for ( int j = 0; j < cell_size; ++j ) {
            // shift de 1 sur l'indexe des noeuds
            connectivity[dim].push_back( ( *nodes_med )[j] - 1 );
        }

        // connectivity_index indique la position des types dans connectivity
        auto szI = connectivity_index[dim].size();
        if ( szI == 0 ) {
            connectivity_index[dim].push_back( 0 );
        }
        connectivity_index[dim].push_back( 1 + connectivity_index[dim].back() + cell_size );

        // Passage entre numerotation globale aster et par dimension medcoupling
        auto sz = corresponding_cells[dim].size();
        corresponding_cells[dim][aster_index - 1] = sz;
    }

    // Tri des groupes d'éléments tojours par dimension
    for ( const auto &group_name : mesh->getGroupsOfCells() ) {
        for ( const auto &cell : mesh->getCells( group_name ) ) {
            int dim = med_to_mc[( *cells_types )[cell]][1];
            groups_c[dim][group_name].push_back( corresponding_cells[dim][cell] );
        }
    }

    // Tri des groupes de noeuds avec shift sur l'indexe des noeuds
    for ( const auto &group_name : mesh->getGroupsOfNodes() ) {
        for ( const auto &aster_node : mesh->getNodes( group_name ) ) {
            groups_n[group_name].push_back( aster_node );
        }
    }

    // Creation des dict Python pour retourner les infos

    PyObject *cells_dict = PyDict_New();

    for ( const auto &iter : connectivity ) {
        auto dim = iter.first;

        PyObject *conn_item = PyTuple_New( 2 );
        PyObject *conn_at_dim = PyTuple_New( connectivity[dim].size() );
        for ( ASTERINTEGER j = 0; j < connectivity[dim].size(); ++j ) {
            PyTuple_SetItem( conn_at_dim, j, PyLong_FromLong( connectivity[dim][j] ) );
        }
        PyTuple_SetItem( conn_item, 0, conn_at_dim );

        PyObject *connI_at_dim = PyTuple_New( connectivity_index[dim].size() );
        for ( ASTERINTEGER j = 0; j < connectivity_index[dim].size(); ++j ) {
            PyTuple_SetItem( connI_at_dim, j, PyLong_FromLong( connectivity_index[dim][j] ) );
        }
        PyTuple_SetItem( conn_item, 1, connI_at_dim );

        PyObject *py_dim = PyLong_FromLong( dim );
        PyDict_SetItem( cells_dict, py_dim, conn_item );

        Py_DECREF( py_dim );
        Py_DECREF( conn_item );
    }

    PyObject *groups_c_dict = PyDict_New();

    for ( const auto &iter : groups_c ) {
        auto dim = iter.first;
        PyObject *groups_c_at_dim = PyDict_New();
        for ( const auto &subiter : groups_c[dim] ) {
            auto gname = subiter.first;
            auto gval = subiter.second;
            PyObject *group_items = PyTuple_New( gval.size() );

            for ( ASTERINTEGER j = 0; j < gval.size(); ++j ) {
                PyTuple_SetItem( group_items, j, PyLong_FromLong( gval[j] ) );
            }

            PyDict_SetItemString( groups_c_at_dim, gname.c_str(), group_items );
            Py_DECREF( group_items );
        }

        PyObject *py_dim = PyLong_FromLong( dim );
        PyDict_SetItem( groups_c_dict, py_dim, groups_c_at_dim );
        Py_DECREF( py_dim );
        Py_DECREF( groups_c_at_dim );
    }

    PyObject *groups_n_dict = PyDict_New();
    for ( const auto &iter : groups_n ) {
        auto gname = iter.first;
        auto gval = iter.second;
        PyObject *group_items = PyTuple_New( gval.size() );
        for ( ASTERINTEGER j = 0; j < gval.size(); ++j ) {
            PyTuple_SetItem( group_items, j, PyLong_FromLong( gval[j] ) );
        }
        PyDict_SetItemString( groups_n_dict, gname.c_str(), group_items );
        Py_DECREF( group_items );
    }

    PyObject *resu_tuple = PyTuple_New( 3 );
    PyTuple_SetItem( resu_tuple, 0, cells_dict );
    PyTuple_SetItem( resu_tuple, 1, groups_c_dict );
    PyTuple_SetItem( resu_tuple, 2, groups_n_dict );
    py::tuple result = py::reinterpret_steal< py::tuple >( resu_tuple );

    return result;
}
