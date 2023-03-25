/**
 * @file ListOfLoads.cxx
 * @brief Implementation de ListOfLoads
 * @author Nicolas Sellenet
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

#include "astercxx.h"

#include "Loads/ListOfLoads.h"

#include "aster_fort_calcul.h"

#include "Supervis/CommandSyntax.h"

#include <typeinfo>

ListOfLoads::ListOfLoads() : ListOfLoads( ModelPtr( NULL ) ) {};

ListOfLoads::ListOfLoads( const std::string &name, const ModelPtr model )
    : DataStructure( name, 19, "L_CHARGES" ),
      _loadInformations( JeveuxVectorLong( getName() + ".INFC" ) ),
      _list( JeveuxVectorChar24( getName() + ".LCHA" ) ),
      _listOfFunctions( JeveuxVectorChar24( getName() + ".FCHA" ) ),
      _isEmpty( true ),
      _model( model ) {};

ListOfLoads::ListOfLoads( const ModelPtr model )
    : ListOfLoads( DataStructureNaming::getNewName( 8 ) + ".LIST_LOAD", model ) {};

ListOfLoads::ListOfLoads( const std::string &name ) : ListOfLoads( name, nullptr ) {};

bool ListOfLoads::checkModelConsistency( const ModelPtr &model ) const {
    if ( _model ) {
        if ( !model ) {
            return false;
        }
        if ( _model != model ) {
            return false;
        }
    }
    return true;
};

bool ListOfLoads::setModel( const ModelPtr &model ) {
    if ( _model ) {
        if ( !this->checkModelConsistency( model ) ) {
            std::string msg =
                "Inconsistent model: " + _model->getName() + " vs " + model->getName();
            AS_ABORT( msg );
        }
    } else
        _model = model;

    return true;
};

int ListOfLoads::getPhysics( void ) const {
    if ( _model )
        return _model->getPhysics();

    return -1;
};

bool ListOfLoads::build( ModelPtr model ) {
    if ( !_isEmpty )
        return true;

    int physic;
    if ( model ) {
        physic = model->getPhysics();
        _model = model;
    } else {
        physic = this->getPhysics();
    }

    ASTERINTEGER iexcit = 1;
    ASTERINTEGER icalc = 0;
    std::string name( getName().c_str() );
    name.resize( 19, ' ' );
    std::string blank( " " );
    blank.resize( 19, ' ' );
    std::string base( JeveuxMemoryTypesNames[(int)Permanent] );

    SyntaxMapContainer dict;
    ListSyntaxMapContainer listeExcit;

    if ( physic == Physics::Mechanics ) {
        CommandSyntax cmdSt( "MECA_STATIQUE" );

        int pos = 0;
        for ( const auto &load : _listOfMechanicalLoadsReal ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = load->getName();
            if ( _listOfMechaFuncReal[pos]->getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfMechaFuncReal[pos]->getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }
        pos = 0;
        for ( const auto &load : _listOfMechanicalLoadsComplex ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = load->getName();
            if ( _listOfMechaFuncComplex[pos]->getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfMechaFuncComplex[pos]->getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }
        pos = 0;
        for ( const auto &load : _listOfMechanicalLoadsFunction ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = load->getName();
            if ( _listOfMechaFuncFunction[pos]->getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfMechaFuncFunction[pos]->getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }
#ifdef ASTER_HAVE_MPI
        pos = 0;
        for ( const auto &load : _listOfParallelMechanicalLoadsReal ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = load->getName();
            if ( _listOfParaMechaFuncReal[pos]->getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfParaMechaFuncReal[pos]->getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }
        pos = 0;
        for ( const auto &load : _listOfParallelMechanicalLoadsFunction ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = load->getName();
            if ( _listOfParaMechaFuncFunction[pos]->getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfParaMechaFuncFunction[pos]->getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }
#endif /* ASTER_HAVE_MPI */

        pos = 0;
        for ( const auto &load : _listOfDirichletBCs ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = load->getName();
            if ( _listOfDiriFun[pos]->getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfDiriFun[pos]->getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }
        dict.container["EXCIT"] = listeExcit;
        cmdSt.define( dict );

        /* Add contact FED */
        std::string ligrel_slav = " ", ligrel_cont = " ";
        if ( _contact.first ) {
            ligrel_slav = _contact.first->getName();
        }
        if ( _contact.second ) {
            ligrel_cont = _contact.second->getName();
        }

        CALLO_NMDOCH_WRAP( name, &iexcit, &icalc, blank, base, ligrel_slav, ligrel_cont );
    } else if ( physic == Physics::Thermal ) {
        CommandSyntax cmdSt( "THER_NON_LINE" );

        int pos = 0;
        for ( const auto &load : _listOfThermalLoadsReal ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = load->getName();
            if ( _listOfTherFuncReal[pos]->getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfTherFuncReal[pos]->getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }

        pos = 0;
        for ( const auto &load : _listOfThermalLoadsFunction ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = load->getName();
            if ( _listOfTherFuncFunction[pos]->getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfTherFuncFunction[pos]->getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }

#ifdef ASTER_HAVE_MPI
        pos = 0;
        for ( const auto &load : _listOfParallelThermalLoadsReal ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = load->getName();
            if ( _listOfParaTherFuncReal[pos]->getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfParaTherFuncReal[pos]->getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }
        pos = 0;
        for ( const auto &load : _listOfParallelThermalLoadsFunction ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = load->getName();
            if ( _listOfParaTherFuncFunction[pos]->getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfParaTherFuncFunction[pos]->getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }
#endif /* ASTER_HAVE_MPI */

        pos = 0;
        for ( const auto &load : _listOfDirichletBCs ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = load->getName();
            if ( _listOfDiriFun[pos]->getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfDiriFun[pos]->getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }

        dict.container["EXCIT"] = listeExcit;
        cmdSt.define( dict );

        CALLO_NTDOCH_WRAP( name, &iexcit, blank, base );
    } else if ( physic == Physics::Acoustic ) {
        CommandSyntax cmdSt( "THER_NON_LINE" );

        int pos = 0;
        for ( const auto &load : _listOfAcousticLoadsComplex ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = load->getName();
            if ( _listOfAcouFuncComplex[pos]->getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfAcouFuncComplex[pos]->getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }

        for ( const auto &load : _listOfDirichletBCs ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = load->getName();
            if ( _listOfDiriFun[pos]->getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfDiriFun[pos]->getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }

        dict.container["EXCIT"] = listeExcit;
        cmdSt.define( dict );

        CALLO_ACDOCH_WRAP( name, &iexcit, blank, base );
    } else {
        AS_ABORT( "Should not be here" );
    }

    _isEmpty = false;
    return true;
};

std::vector< FiniteElementDescriptorPtr > ListOfLoads::getFiniteElementDescriptors() const {
    std::vector< FiniteElementDescriptorPtr > FEDesc;

    const int physic = this->getPhysics();

    if ( physic == Physics::Mechanics ) {
        for ( const auto &load : _listOfMechanicalLoadsReal ) {
            FEDesc.push_back( load->getFiniteElementDescriptor() );
        }
        for ( const auto &load : _listOfMechanicalLoadsComplex ) {
            FEDesc.push_back( load->getFiniteElementDescriptor() );
        }
        for ( const auto &load : _listOfMechanicalLoadsFunction ) {
            FEDesc.push_back( load->getFiniteElementDescriptor() );
        }
#ifdef ASTER_HAVE_MPI
        for ( const auto &load : _listOfParallelMechanicalLoadsReal ) {
            FEDesc.push_back( load->getFiniteElementDescriptor() );
        }
        for ( const auto &load : _listOfParallelMechanicalLoadsFunction ) {
            FEDesc.push_back( load->getFiniteElementDescriptor() );
        }
#endif /* ASTER_HAVE_MPI */
    } else if ( physic == Physics::Thermal ) {
        for ( const auto &load : _listOfThermalLoadsReal ) {
            FEDesc.push_back( load->getFiniteElementDescriptor() );
        }
        for ( const auto &load : _listOfThermalLoadsFunction ) {
            FEDesc.push_back( load->getFiniteElementDescriptor() );
        }
    } else if ( physic == Physics::Acoustic ) {
        for ( const auto &load : _listOfAcousticLoadsComplex ) {
            FEDesc.push_back( load->getFiniteElementDescriptor() );
        }
#ifdef ASTER_HAVE_MPI
        for ( const auto &load : _listOfParallelThermalLoadsReal ) {
            FEDesc.push_back( load->getFiniteElementDescriptor() );
        }
        for ( const auto &load : _listOfParallelThermalLoadsFunction ) {
            FEDesc.push_back( load->getFiniteElementDescriptor() );
        }
#endif /* ASTER_HAVE_MPI */
    }

    return FEDesc;
};
