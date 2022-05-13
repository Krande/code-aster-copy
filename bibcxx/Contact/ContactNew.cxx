/**
 * @file ContactNew.cxx
 * @brief Implementation de Contact
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

#include "Contact/ContactNew.h"

#include "aster_fort_ds.h"
#include "aster_fort_jeveux.h"
#include "aster_fort_mesh.h"
#include "aster_fort_utils.h"

#include "Messages/Messages.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Utilities/Tools.h"

using VectorLongIter = VectorLong::iterator;

ContactNew::ContactNew( const std::string name, const ModelPtr model )
    : DataStructure( name, 8, "CHAR_CONT" ),
      _model( model ),
      _FEDesc( std::make_shared< FiniteElementDescriptor >( getName() + ".CHME.LIGRE",
                                                            _model->getMesh() ) ),
      _verbosity( 1 ),
      _friction( false ),
      _smoothing( false ) {
    // model has to be mechanics
    if ( !_model->isMechanical() )
        UTMESS( "F", "CONTACT1_2" );
};

bool ContactNew::build() {
    auto mesh = getMesh();
    ASTERINTEGER nb_zones = getNumberOfContactZones();
    ASTERINTEGER nb_slave_cells = 0;

    // sdcont_defi.CONTACT.MAILCO/NOEUCO/ssnoco
    std::vector< std::pair< VectorLong, VectorLong > > mailco;
    std::vector< std::pair< VectorLong, VectorLong > > noeuco;

    for ( long i = 0; i < nb_zones; i++ ) {
        auto zone_i = getContactZone( i );
        std::string slave = zone_i->getSlaveGroupOfCells();
        std::string master = zone_i->getMasterGroupOfCells();
        VectorString sans_grs = zone_i->getExcludedSlaveGroupOfCells();

        // read slave/master nodes/cell : localNumbering, same_rank
        auto l_slave_nodes = mesh->getNodesFromCells( slave, false, true );
        auto l_slave_cells = mesh->getCells( slave );
        auto l_master_nodes = mesh->getNodesFromCells( master, false, true );
        auto l_master_cells = mesh->getCells( master );

        // read hte nodes by SANS_GROUP_MA : several groups
        VectorLong l_sans_nodes;
        for ( auto &name : sans_grs ) {
            VectorLongIter it = l_sans_nodes.end();
            VectorLong sans_gr_i = mesh->getNodesFromCells( name, false, true );
            l_sans_nodes.insert( it, sans_gr_i.begin(), sans_gr_i.end() );
        }

        // check between SANS_GROUP_MA and slave nodes : save in commonNodes
        VectorLong commonNodes;
        if ( mesh->isParallel() ) {
#ifdef ASTER_HAVE_MPI
            VectorLong lg_sans_nodes;
            AsterMPI::all_gather( l_sans_nodes, lg_sans_nodes );
            commonNodes = set_intersection( lg_sans_nodes, l_slave_nodes );
#endif
        } else {
            commonNodes = set_intersection( l_sans_nodes, l_slave_nodes );
        }
        // save info
        nb_slave_cells = nb_slave_cells + l_slave_cells.size();
        mailco.push_back( std::make_pair( l_master_cells, l_slave_cells ) );
        noeuco.push_back( std::make_pair( l_master_nodes, l_slave_nodes ) );
    }

    // check the common slave nodes between les zones
    VectorLong doublNodes;
    for ( auto it = noeuco.begin(); it != noeuco.end(); ++it ) {
        VectorLong l_a = it->second;
        for ( auto itb = std::next( it, 1 ); itb != noeuco.end(); ++itb ) {
            VectorLong l_b = itb->second;
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
                nb_doublNodes = 0;
                AsterMPI::all_reduce( nb_doublNodes_lc, nb_doublNodes, MPI_MIN );
            }
#endif
            if ( nb_doublNodes > 0 ) {
                for(auto& node: doublNodes)
                    std::cout << node << std::endl;
                UTMESS( "F", "CONTACT1_3" );
            }
        }
    }

    // create slave elements in model : routine mmprel
    std::string ligret = ljust( "&&OP0030.LIGRET", 19, ' ' );
    std::string phenom = ljust( "MECANIQUE", 16, ' ' );
    std::string modeli;

    // calico
    ASTERINTEGER nb_dim = 0;
    ASTERINTEGER nb_dim_ = _model->getGeometricDimension();
    if ( nb_dim_ > 3 ) { // general ? dans model ?
        UTMESS( "A", "CONTACT_84" );
        if ( nb_dim_ == 1003 ) {
            nb_dim = 3;
        } else if ( nb_dim_ == 1002 ) {
            nb_dim = 2;
        } else if ( nb_dim_ == 23 ) {
            nb_dim = 2;
        } else {
            UTMESS( "F", "CONTACT1_4" );
        }
    } else {
        nb_dim = nb_dim_;
    }

    CALL_JEMARQ();
    ASTERINTEGER i_zone = 0;
    std::string jeveuxname = ljust( "&&MMPREL.LISTE_MAILLES", 24, ' ' );
    for ( auto &[mastercells, slavecells] : mailco ) {
        bool hasFr = getContactZone( i_zone )->getFrictionParameter()->hasFriction();
        if ( nb_dim == 2 ) {
            if ( hasFr ) {
                modeli = "FRIC_SL_2D";
            } else {
                modeli = "CONT_SL_2D";
            }
        } else if ( nb_dim == 3 ) {
            if ( hasFr ) {
                modeli = "FRIC_SL_3D";
            } else {
                modeli = "CONT_SL_3D";
            }
        }
        modeli = ljust( modeli, 16, ' ' );

        // AFFECTATION DE L'OBJET DE TYPE LIGRET ET DE NOM LIGREZ jeveuxname
        ASTERINTEGER slave_cells_i = slavecells.size();
        if ( slave_cells_i > 0 ) {
            for (auto &cell : slavecells)
                cell += 1;
            JeveuxVectorLong list_elem = JeveuxVectorLong( jeveuxname, slavecells );
            CALL_AJELLT( ligret.c_str(), mesh->getName().c_str(), &slave_cells_i,
                         jeveuxname.c_str(), " ", phenom.c_str(), modeli.c_str(), 0, " " );
        }
        i_zone++;
    }
    CALL_JEDETR( jeveuxname.c_str() );
    CALL_JEDEMA();

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
    }

    return true;
}
