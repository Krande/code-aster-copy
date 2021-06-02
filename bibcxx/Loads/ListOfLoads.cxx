/**
 * @file ListOfLoads.cxx
 * @brief Implementation de ListOfLoads
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

#include <stdexcept>
#include <typeinfo>
#include "astercxx.h"

#include "aster_fort_calcul.h"
#include "Loads/ListOfLoads.h"
#include "Supervis/CommandSyntax.h"

ListOfLoadsClass::ListOfLoadsClass( const JeveuxMemory memType )
    : DataStructure( DataStructureNaming::getNewName( memType, 8 ) + ".LIST_LOAD", 19, "L_CHARGES",
                     memType ),
      _loadInformations( JeveuxVectorLong( getName() + ".INFC" ) ),
      _list( JeveuxVectorChar24( getName() + ".LCHA" ) ),
      _listOfFunctions( JeveuxVectorChar24( getName() + ".FCHA" ) ),
      _isEmpty( true ), _model( nullptr ){};


bool ListOfLoadsClass::checkModelConsistency( const ModelPtr& model ) const
{
    if( _model ){
        if( _model->getName() != model->getName())
            return false;
    }
    return true;
};

bool ListOfLoadsClass::setModel( const ModelPtr& model )
{
    if( _model ){
        if( !this->checkModelConsistency( model ) )
            throw std::runtime_error("Inconsistent model");
    }
    else
        _model = model;

    return true;
};

int ListOfLoadsClass::getPhysics( void ) const
{
    if( _model )
        return _model->getPhysics();

    return -1;
};

bool ListOfLoadsClass::build( ModelPtr model ) {
    if ( !_isEmpty )
        return true;

    int physic;
    if( model )
        physic = model->getPhysics();
    else
        physic = this->getPhysics();

    ASTERINTEGER iexcit = 1;
    std::string name( getName().c_str() );
    name.resize( 19, ' ' );
    std::string blank( " " );
    blank.resize( 19, ' ' );
    std::string base( JeveuxMemoryTypesNames[(int)getMemoryType()] );

    SyntaxMapContainer dict;
    ListSyntaxMapContainer listeExcit;

    if( physic == Physics::Mechanics)
    {
        CommandSyntax cmdSt( "MECA_STATIQUE" );
        int pos = 0;
        for ( const auto &curIter : _listOfMechanicalLoadsReal ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = curIter->getName();
            if ( _listOfMechaFuncReal[pos].getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfMechaFuncReal[pos].getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }
        pos = 0;
        for ( const auto &curIter : _listOfMechanicalLoadsFunction ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = curIter->getName();
            if ( _listOfMechaFuncFunction[pos].getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfMechaFuncFunction[pos].getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }
#ifdef ASTER_HAVE_MPI
        pos = 0;
        for ( const auto &curIter : _listOfParallelMechanicalLoadsReal ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = curIter->getName();
            if ( _listOfParaMechaFuncReal[pos].getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfParaMechaFuncReal[pos].getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }
        pos = 0;
        for ( const auto &curIter : _listOfParallelMechanicalLoadsFunction ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = curIter->getName();
            if ( _listOfParaMechaFuncFunction[pos].getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfParaMechaFuncFunction[pos].getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }
#endif /* ASTER_HAVE_MPI */

        pos = 0;
        for ( const auto &curIter : _listOfDirichletBCs ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = curIter->getName();
            if ( _listOfDiriFun[pos].getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfDiriFun[pos].getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }
        dict.container["EXCIT"] = listeExcit;
        cmdSt.define( dict );

        CALLO_NMDOCH_WRAP( name, &iexcit, blank, base );
    }
    else if( physic == Physics::Thermal)
    {
        CommandSyntax cmdSt( "THER_NON_LINE" );
        int pos = 0;
        for ( const auto &curIter : _listOfThermalLoads ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = curIter->getName();
            if ( _listOfTherFun[pos].getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfTherFun[pos].getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }

        pos = 0;
        for ( const auto &curIter : _listOfDirichletBCs ) {
            SyntaxMapContainer dict2;
            dict2.container["CHARGE"] = curIter->getName();
            if ( _listOfDiriFun[pos].getName() != emptyRealFunction->getName() )
                dict2.container["FONC_MULT"] = _listOfDiriFun[pos].getName();
            ++pos;
            listeExcit.push_back( dict2 );
        }

        dict.container["EXCIT"] = listeExcit;
        cmdSt.define( dict );

        CALLO_NTDOCH_WRAP( name, &iexcit, blank, base );
    }
    else if( physic == Physics::Acoustic)
    {
        throw std::runtime_error("Not Implemented AcousticLoad");
    }
    else{
        throw std::runtime_error("Should not be here");
    }

    _isEmpty = false;
    return true;
};


std::vector< FiniteElementDescriptorPtr > ListOfLoadsClass::getFiniteElementDescriptors() const
{
    std::vector< FiniteElementDescriptorPtr > FEDesc;

    const int physic = this->getPhysics();

    if( physic == Physics::Mechanics)
    {
        for ( const auto &curIter : _listOfMechanicalLoadsReal ) {
            FEDesc.push_back( curIter->getFiniteElementDescriptor() );
        }
        for ( const auto &curIter : _listOfMechanicalLoadsFunction ) {
            FEDesc.push_back( curIter->getFiniteElementDescriptor() );
        }
#ifdef ASTER_HAVE_MPI
        for ( const auto &curIter : _listOfParallelMechanicalLoadsReal ) {
            FEDesc.push_back( curIter->getFiniteElementDescriptor() );
        }
        for ( const auto &curIter : _listOfParallelMechanicalLoadsFunction ) {
            FEDesc.push_back( curIter->getFiniteElementDescriptor() );
        }
#endif /* ASTER_HAVE_MPI */
    }
    else if( physic == Physics::Thermal)
    {
        for ( const auto &curIter : _listOfThermalLoads ) {
            FEDesc.push_back( curIter->getFiniteElementDescriptor() );
        }
    }
    else if( physic == Physics::Acoustic)
    {
        for ( const auto &curIter : _listOfAcousticLoads ) {
            FEDesc.push_back( curIter->getFiniteElementDescriptor() );
        }
    }

    return FEDesc;
};
