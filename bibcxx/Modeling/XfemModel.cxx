/**
 * @file XfemModel.cxx
 * @brief Implementation de Model
 * @author
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

#include "Modeling/XfemModel.h"

#include "DataFields/FieldOnNodes.h"
#include "DataFields/FieldOnCells.h"


XfemModel::SubElementTopology::SubElementTopology( const std::string name):
                            _name( name + ".TOPOSE" ),
                            pin( new FieldOnCellsReal( getName() + ".PIN" ) ),
                            cns( new FieldOnCellsLong( getName() + ".CNS" ) ),
                            hea( new FieldOnCellsLong( getName() + ".HEA" ) ),
                            lon( new FieldOnCellsLong( getName() + ".LON" ) ),
                            pai( new FieldOnCellsReal( getName() + ".PAI" ) ),
                            pmi( new FieldOnCellsReal( getName() + ".PMI" ) ){};

XfemModel::FacetTopology::FacetTopology( const std::string name):
                            _name( name + ".TOPOFAC" ),
                            _intersection_pt( new FieldOnCellsReal( getName() + ".PI" ) ),
                            _intersection_edge( new FieldOnCellsReal( getName() + ".AI" ) ),
                            _connectivity( new FieldOnCellsLong( getName() + ".CF" ) ),
                            _length( new FieldOnCellsLong( getName() + ".LO" ) ),
                            _base( new FieldOnCellsReal( getName() + ".BA" ) ),
                            _heaviside( new FieldOnCellsLong( getName() + ".HE" ) ),
                            _intersection_pt2( new FieldOnCellsReal( getName() + ".OE" ) ){};

XfemModel::NodalTopology::NodalTopology( const std::string name):
                            _name( name + ".TOPONO" ),
                            hno( new FieldOnCellsLong( getName() + ".HNO" ) ),
                            hfa( new FieldOnCellsLong( getName() + ".HFA" ) ),
                            hse( new FieldOnCellsLong( getName() + ".HSE" ) ){};

XfemModel::XfemModel( const std::string name ):
                                _topose( SubElementTopology( name ) ),
                                _topofac( FacetTopology( name ) ),
                                _topono( NodalTopology( name ) ),
                                _normal_levelset( new FieldOnCellsReal( name + ".LNNO" ) ),
                                _tangent_levelset( new FieldOnCellsReal( name + ".LTNO") ),
                                _local_basis( new FieldOnCellsReal( name + ".BASLOC" ) ),
                                _nodal_status( new FieldOnCellsLong( name + ".STNO" ) ),
                                _crack_nodes( new FieldOnCellsLong( name + ".FISSNO" ) ),
                                _crack_conn( new FieldOnCellsLong( name + ".FISSCO" ) ),
                                _heaviside( new FieldOnCellsLong( name + ".HEAVNO" ) ),
                                _cracked_cells( new FieldOnCellsLong( name + ".XMAFIS" ) ),
                                _xfem_nodes( new FieldOnNodesLong( name + ".NOXFEM" ) ),
                                _contact( JeveuxVectorLong( name + ".XFEM_CONT" ) ),
                                _crack_number( JeveuxVectorLong( name + ".NFIS" ) ),
                                _crack_names( JeveuxVectorChar8( name + ".FISS" ) ),
                                _pre_cond( JeveuxVectorChar8( name + ".PRE_COND" ) ),
                                _thermic( JeveuxVectorChar8( name + ".MODELE_THER" ) ){
                                    _listfields.insert( { "PINTTO", _topose.pin } );
                                    _listfields.insert( { "CNSETO", _topose.cns } );
                                    _listfields.insert( { "HEAVTO", _topose.hea } );
                                    _listfields.insert( { "LONCHA", _topose.lon } );
                                    _listfields.insert( { "PMILT", _topose.pmi });
                                    _listfields.insert( { "HEAVTO", _topose.hea } );  
                                    _listfields.insert( { "BASLOC", _local_basis } );
                                    _listfields.insert( { "LSN", _normal_levelset } );
                                    _listfields.insert( { "LST", _tangent_levelset } );
                                    _listfields.insert( { "STANO", _nodal_status } );
                                    _listfields.insert( { "FISSNO", _crack_nodes } );  
                                    _listfields.insert( { "HEAVNO", _topono.hno } );
                                    _listfields.insert( { "HEAVFA", _topono.hfa } );
};
