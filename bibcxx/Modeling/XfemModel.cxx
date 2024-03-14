/**
 * @file XfemModel.cxx
 * @brief Implementation de Model
 * @author
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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

#include "Modeling/XfemModel.h"

#include "DataFields/FieldOnCells.h"
#include "DataFields/FieldOnNodes.h"

XfemModel::SubElementTopology::SubElementTopology( const std::string name )
    : _name( name + ".TOPOSE" ),
      pin( std::make_shared< FieldOnCellsReal >( getName() + ".PIN" ) ),
      cns( std::make_shared< FieldOnCellsLong >( getName() + ".CNS" ) ),
      hea( std::make_shared< FieldOnCellsLong >( getName() + ".HEA" ) ),
      lon( std::make_shared< FieldOnCellsLong >( getName() + ".LON" ) ),
      pai( std::make_shared< FieldOnCellsReal >( getName() + ".PAI" ) ),
      pmi( std::make_shared< FieldOnCellsReal >( getName() + ".PMI" ) ) {};

XfemModel::FacetTopology::FacetTopology( const std::string name )
    : _name( name + ".TOPOFAC" ),
      intersection_pt( std::make_shared< FieldOnCellsReal >( getName() + ".PI" ) ),
      intersection_edge( std::make_shared< FieldOnCellsReal >( getName() + ".AI" ) ),
      connectivity( std::make_shared< FieldOnCellsLong >( getName() + ".CF" ) ),
      length( std::make_shared< FieldOnCellsLong >( getName() + ".LO" ) ),
      base( std::make_shared< FieldOnCellsReal >( getName() + ".BA" ) ),
      _heaviside( std::make_shared< FieldOnCellsLong >( getName() + ".HE" ) ),
      intersection_pt2( std::make_shared< FieldOnCellsReal >( getName() + ".OE" ) ) {};

XfemModel::NodalTopology::NodalTopology( const std::string name )
    : _name( name + ".TOPONO" ),
      hno( std::make_shared< FieldOnCellsLong >( getName() + ".HNO" ) ),
      hfa( std::make_shared< FieldOnCellsLong >( getName() + ".HFA" ) ),
      hse( std::make_shared< FieldOnCellsLong >( getName() + ".HSE" ) ) {};

XfemModel::XfemModel( const std::string name )
    : _topose( SubElementTopology( name ) ),
      _topofac( FacetTopology( name ) ),
      _topono( NodalTopology( name ) ),
      _normal_levelset( std::make_shared< FieldOnCellsReal >( name + ".LNNO" ) ),
      _tangent_levelset( std::make_shared< FieldOnCellsReal >( name + ".LTNO" ) ),
      _local_basis( std::make_shared< FieldOnCellsReal >( name + ".BASLOC" ) ),
      _nodal_status( std::make_shared< FieldOnCellsLong >( name + ".STNO" ) ),
      _crack_nodes( std::make_shared< FieldOnCellsLong >( name + ".FISSNO" ) ),
      _crack_conn( std::make_shared< FieldOnCellsLong >( name + ".FISSCO" ) ),
      _heaviside( std::make_shared< FieldOnCellsLong >( name + ".HEAVNO" ) ),
      _cracked_cells( std::make_shared< FieldOnCellsLong >( name + ".XMAFIS" ) ),
      _xfem_nodes( std::make_shared< FieldOnNodesLong >( name + ".NOXFEM" ) ),
      _contact( JeveuxVectorLong( name + ".XFEM_CONT" ) ),
      _crack_number( JeveuxVectorLong( name + ".NFIS" ) ),
      _crack_names( JeveuxVectorChar8( name + ".FISS" ) ),
      _pre_cond( JeveuxVectorChar8( name + ".PRE_COND" ) ),
      _thermic( JeveuxVectorChar8( name + ".MODELE_THER" ) ) {
    _xfem_nodes->setDescription( std::make_shared< EquationNumbering >( name + ".NOXF.NUMEQ" ) );
    _listfields.insert( { "PINTTO", _topose.pin } );
    _listfields.insert( { "CNSETO", _topose.cns } );
    _listfields.insert( { "HEAVTO", _topose.hea } );
    _listfields.insert( { "LONCHA", _topose.lon } );
    _listfields.insert( { "PMILT", _topose.pmi } );
    _listfields.insert( { "HEAVTO", _topose.hea } );
    _listfields.insert( { "BASLOC", _local_basis } );
    _listfields.insert( { "LSN", _normal_levelset } );
    _listfields.insert( { "LST", _tangent_levelset } );
    _listfields.insert( { "STANO", _nodal_status } );
    _listfields.insert( { "FISSNO", _crack_nodes } );
    _listfields.insert( { "HEAVNO", _topono.hno } );
    _listfields.insert( { "HEAVSE", _topono.hse } );
    _listfields.insert( { "HEAVFA", _topono.hfa } );
    _listfields.insert( { "AINTER", _topofac.intersection_edge } );
    _listfields.insert( { "PINTER", _topofac.intersection_pt2 } );
    _listfields.insert( { "CFACE", _topofac.connectivity } );
    _listfields.insert( { "LONGCO", _topofac.length } );
    _listfields.insert( { "BASECO", _topofac.base } );
};

ASTERINTEGER XfemModel::getContact() const {
    _contact->updateValuePointer();
    return ( *_contact )[0];
};
