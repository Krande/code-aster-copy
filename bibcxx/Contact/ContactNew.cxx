/**
 * @file ContactNew.cxx
 * @brief Implementation de Contact
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

#include "Contact/ContactNew.h"

#include "aster_fort_ds.h"
#include "aster_fort_jeveux.h"
#include "aster_fort_mesh.h"
#include "aster_fort_utils.h"

#include "Meshes/ConnectionMesh.h"
#include "Messages/Messages.h"
#include "Modeling/Model.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Utilities/Tools.h"

using VectorLongIter = VectorLong::iterator;

ContactNew::ContactNew( const std::string name, const ModelPtr model, const std::string type )
    : DSWithCppPickling( name, 8, type ), _model( model ), _FEDesc( nullptr ), _verbosity( 1 ) {};

/** @brief restricted constructor (Set) and method (Get) to support pickling */
ContactNew::ContactNew( const py::tuple &tup )
    : ContactNew( tup[0].cast< std::string >(), tup[1].cast< ModelPtr >() ) {
    int i = 1;
    _FEDesc = tup[++i].cast< FiniteElementDescriptorPtr >();
    _zones = tup[++i].cast< std::vector< ContactZonePtr > >();
    _verbosity = tup[++i].cast< ASTERINTEGER >();
    // _FEDesc =
    //     std::make_shared< FiniteElementDescriptor >( getName() + ".CONT.LIGRE", _model->getMesh()
    //     );
}
py::tuple ContactNew::_getState() const {
    return py::make_tuple( getName(), _model, _FEDesc, _zones, _verbosity );
}

FrictionNew::FrictionNew( const py::tuple &tup )
    : FrictionNew( tup[0].cast< std::string >(), tup[1].cast< ModelPtr >() ) {
    int i = 1;
    _FEDesc = tup[++i].cast< FiniteElementDescriptorPtr >();
    _zones = tup[++i].cast< std::vector< ContactZonePtr > >();
    _verbosity = tup[++i].cast< ASTERINTEGER >();
}

py::tuple FrictionNew::_getState() const {
    return py::make_tuple( getName(), _model, _FEDesc, _zones, _verbosity );
}

void ContactNew::appendContactZone( const ContactZonePtr zone ) {
    _zones.push_back( zone );
    _zones.back()->setVerbosity( getVerbosity() );
}

ASTERINTEGER ContactNew::getSpaceDime() const {

    ASTERINTEGER spaceDime = 0;
    ASTERINTEGER spaceDime_ = getModel()->getGeometricDimension();
    if ( spaceDime_ > 3 ) {
        UTMESS( "A", "CONTACT_84" );
        if ( spaceDime_ == 1003 ) {
            spaceDime = 3;
        } else if ( spaceDime_ == 1002 ) {
            spaceDime = 2;
        } else if ( spaceDime_ == 23 ) {
            spaceDime = 2;
        } else {
            UTMESS( "F", "CONTACT1_4" );
        }
    } else {
        spaceDime = spaceDime_;
    }
    return spaceDime;
}

bool ContactNew::build() {
    CALL_JEMARQ();
    // model has to be mechanics
    if ( !_model->isMechanical() )
        UTMESS( "F", "CONTACT1_2" );
    _FEDesc =
        std::make_shared< FiniteElementDescriptor >( getName() + ".CONT.LIGRE", _model->getMesh() );

    auto mesh = getMesh();

    // Space dimension
    ASTERINTEGER spaceDime = this->getSpaceDime();

    // Define name of catalogue
    std::map< std::tuple< ASTERINTEGER, ContactAlgo, bool >, std::string > cata;

    cata[std::make_tuple( 2, ContactAlgo::Lagrangian, false )] = "CONT_LAG_SL_2D";
    cata[std::make_tuple( 3, ContactAlgo::Lagrangian, false )] = "CONT_LAG_SL_3D";
    cata[std::make_tuple( 2, ContactAlgo::Lagrangian, true )] = "FRIC_LAG_SL_2D";
    cata[std::make_tuple( 3, ContactAlgo::Lagrangian, true )] = "FRIC_LAG_SL_3D";

    // same cata for contact and friction since no Lagrange
    cata[std::make_tuple( 2, ContactAlgo::Nitsche, false )] = "CONT_NIT_SL_2D";
    cata[std::make_tuple( 3, ContactAlgo::Nitsche, false )] = "CONT_NIT_SL_3D";
    cata[std::make_tuple( 2, ContactAlgo::Nitsche, true )] = "CONT_NIT_SL_2D";
    cata[std::make_tuple( 3, ContactAlgo::Nitsche, true )] = "CONT_NIT_SL_3D";

    // Same elements for contact and friction
    cata[std::make_tuple( 2, ContactAlgo::Penalization, false )] = "CONT_PENA_SL_2D";
    cata[std::make_tuple( 3, ContactAlgo::Penalization, false )] = "CONT_PENA_SL_3D";
    cata[std::make_tuple( 2, ContactAlgo::Penalization, true )] = "CONT_PENA_SL_2D";
    cata[std::make_tuple( 3, ContactAlgo::Penalization, true )] = "CONT_PENA_SL_3D";

    ASTERINTEGER nb_slave_cells = 0;

    // List of all contact cells
    std::vector< std::pair< VectorLong, std::string > > mailco;

    // List of all contact nodes
    std::vector< VectorLong > noeuco;

    for ( auto &zone_i : _zones ) {
        zone_i->build( _model );
        // read slave nodes/cell : localNumbering, same_rank
        auto l_slave_nodes = zone_i->getSlaveNodes();
        auto l_slave_cells = zone_i->getSlaveCells();

        // save info
        nb_slave_cells += l_slave_cells.size();
        mailco.push_back( std::make_pair(
            l_slave_cells,
            cata[std::make_tuple( spaceDime, zone_i->getContactParameter()->getAlgorithm(),
                                  zone_i->hasFriction() )] ) );
        noeuco.push_back( l_slave_nodes );
    }

    // Check the common slave nodes between zones (or cells ?)
    VectorLong doublNodes;
    for ( auto it = noeuco.begin(); it != noeuco.end(); ++it ) {
        VectorLong l_a = *it;
        for ( auto itb = std::next( it, 1 ); itb != noeuco.end(); ++itb ) {
            VectorLong l_b = *itb;
            if ( mesh->isParallel() ) {
#ifdef ASTER_HAVE_MPI
                VectorLong lg_a;
                AsterMPI::all_gather( l_a, lg_a );
                doublNodes = set_intersection( lg_a, l_b );
#endif
            } else {
                doublNodes = set_intersection( l_a, l_b );
            }
            ASTERINTEGER nb_doublNodes = doublNodes.size();
// share error
#ifdef ASTER_HAVE_MPI
            if ( mesh->isParallel() ) {
                ASTERINTEGER nb_doublNodes_lc = nb_doublNodes;
                nb_doublNodes = AsterMPI::max( nb_doublNodes_lc );
            }
#endif
            if ( nb_doublNodes > 0 ) {
                UTMESS( "F", "CONTACT1_1" );
            }
        }
    }

    // Create slave elements in model : routine mmprel
    std::string ligret = ljust( "&&OP0030.LIGRET", 19, ' ' );
    std::string phenom = ljust( "MECANIQUE", 16, ' ' );
    std::string modeli;

    std::string jeveuxname = ljust( "&&MMPREL.LISTE_MAILLES", 24, ' ' );
    for ( auto &[slavecells, modeli] : mailco ) {
        modeli = ljust( modeli, 16, ' ' );

        // AFFECTATION DE L'OBJET DE TYPE LIGRET ET DE NOM LIGREZ jeveuxname
        ASTERINTEGER slave_cells_i = slavecells.size();
        if ( slave_cells_i > 0 ) {
            std::transform( slavecells.begin(), slavecells.end(), slavecells.begin(),
                            []( ASTERINTEGER &i ) -> ASTERINTEGER { return ++i; } );

            JeveuxVectorLong list_elem = JeveuxVectorLong( jeveuxname, slavecells );
            CALL_AJELLT( ligret.c_str(), mesh->getName().c_str(), &slave_cells_i,
                         list_elem->getDataPtr(), phenom.c_str(), modeli.c_str() );
        }
    }

    // hpc : only when the proc has slave cells, so with ligret
    if ( nb_slave_cells != 0 ) {
        std::string ligrel = ljust( "&&OP0030.LIGREL", 19, ' ' );
        std::string base = "V";
        CALL_LGTLGR( base.c_str(), ligret.c_str(), ligrel.c_str() );

        std::string type = "LIGRET";
        CALLO_DETRSD( type, ligret );

        base = "G";
        type = "LIGREL";
        CALLO_COPISD( type, base, ligrel, _FEDesc->getName() );
        CALLO_DETRSD( type, ligrel );

        CALLO_ADALIG_WRAP( _FEDesc->getName() );
        CALLO_CORMGI( base, _FEDesc->getName() );

        std::string _name = _FEDesc->getName() + ".LGRF";
        std::string param = "DOCU";
        std::string value = "MECA";
        CALLO_JEECRA_STRING_WRAP( _name, param, value );
        bool l_calc_rigi = false;
        CALLO_INITEL( _FEDesc->getName(), (ASTERLOGICAL *)&l_calc_rigi );

        _FEDesc->build();
    }

    CALL_JEDEMA();

    return true;
}

VectorLong ContactNew::getSlaveNodes() const {

    SetLong nodes;

    for ( auto &zone : _zones ) {
        auto l_slave_nodes = zone->getSlaveNodes();

        nodes.insert( l_slave_nodes.begin(), l_slave_nodes.end() );
    }

    return VectorLong( nodes.begin(), nodes.end() );
};

VectorLong ContactNew::getSlaveCells() const {

    SetLong cells;

    for ( auto &zone : _zones ) {
        auto l_slave_cells = zone->getSlaveCells();

        cells.insert( l_slave_cells.begin(), l_slave_cells.end() );
    }
    return VectorLong( cells.begin(), cells.end() );
};

void ContactNew::enableFriction( const bool &friction ) {
    for ( auto &zone : _zones ) {
        zone->enableFriction( friction );
    }
};

bool ContactNew::hasFriction() const {
    for ( auto &zone : _zones ) {
        if ( zone->hasFriction() ) {
            return true;
        }
    }
    return false;
};

void ContactNew::enableSmoothing( const bool &smoothing ) {
    for ( auto &zone : _zones ) {
        zone->enableSmoothing( smoothing );
    }
};

bool ContactNew::hasSmoothing() const {
    for ( auto &zone : _zones ) {
        if ( zone->hasSmoothing() ) {
            return true;
        }
    }
    return false;
};

void ContactNew::setVerbosity( const ASTERINTEGER &level ) {
    _verbosity = level;
    for ( auto &zone : _zones ) {
        zone->setVerbosity( level );
    }
}

VectorLong ContactNew::getNumberOfIntersectionPoints() const {
    VectorLong returnValue;

    return returnValue;
}

VectorLong ContactNew::getNumberOfIntersectionPoints( const ASTERINTEGER &indexZone ) const {
    VectorLong returnValue;

    return returnValue;
}
