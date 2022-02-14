/**
 * @file ContactPairing.cxx
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

#include "Contact/ContactPairing.h"




ContactPairing::ContactPairing( const std::string name, const std::vector< ContactZonePtr > zones,
     const BaseMeshPtr mesh): DataStructure( name, 8, "PAIRING_SD" ), _zones(zones), _mesh(mesh)
    {
      _newCoordinates = boost::make_shared<MeshCoordinatesField>(*(_mesh->getCoordinates()));

     // be sure that zones is not empty and get size of zones and resize
     if(_zones.empty()) throw std::runtime_error(" ContactZone object is empty ");
     int size_zones = zones.size();

      // resize pairing quantities
      _nbPairs.resize(size_zones);
      _listOfPairs.resize(size_zones);
      _nbIntersectionPoints.resize(size_zones);
      _slaveIntersectionPoints.resize(size_zones);
      _masterIntersectionPoints.resize(size_zones);
      _gaussPoints.resize(size_zones);

    };


 void ContactPairing::updateCoordinates(FieldOnNodesRealPtr& disp){

      std::string base("V"), cumul("CUMU");
      ASTERDOUBLE alpha = 1.;

      MeshCoordinatesFieldPtr oldCoord =  _mesh->getCoordinates();


      CALLO_VTGPLD(cumul, &alpha, oldCoord->getName(), disp->getName(), 
                                        base, _newCoordinates->getName());

      try{
          _newCoordinates->updateValuePointers();
      }catch(const std::exception& e){
           std::cout << e.what() << std::endl;
      }
      
 }


ASTERBOOL ContactPairing::computePairingQuantities( ASTERINTEGER i ){

     if(i < 0 || i >= _zones.size()) {
          throw std::out_of_range( "The zone index should be between 0  and " 
                              + std::to_string(_zones.size() - 1) );
     }


     // get and define some input parameters
     VectorLong&  eleMaster    =  _zones[i]->getMasterCells();
     VectorLong&  NodesMaster  =  _zones[i]->getMasterNodes();
     VectorLong&  eleSlave     =  _zones[i]->getSlaveCells();
     ASTERINTEGER nbCellMaster =  eleMaster.size();
     ASTERINTEGER nbNodeMaster =  NodesMaster.size();
     ASTERINTEGER nbCellSlave  =  eleSlave.size();
     std::string pair_method;

     // get pairing method
     ContactVariant variant = _zones[i]->getContactParameter()->getVariant();
     if( variant == ContactVariant::Robust ){
          pair_method = ljust( "ROBUSTE", 24, ' ' );
     }else if( variant == ContactVariant::Rapide ){
          pair_method = ljust( "RAPIDE", 24, ' ' );
     }else {
          pair_method = ljust( "ROBUSTE", 24, ' ' );
     }

     // to remove later
     if( variant != ContactVariant::Robust) std::cout << "WARNING: ONLY "
                                   "ROBUST VARIANT IS AVAILBALE" << std::endl;
     pair_method = ljust( "ROBUSTE", 24, ' ' );

     // tolerence
     ASTERDOUBLE pair_tole = 0.001;
     
     // set pairs numbers to 0
     ASTERINTEGER nb_pairs = 0;

     // output paramaters as C pointers
     ASTERINTEGER *pairs = NULL;
     ASTERINTEGER *nbInterPoints = NULL;
     ASTERDOUBLE  *InterSlavePoints = NULL;
     ASTERDOUBLE  *InterMasterPoints = NULL;
     ASTERDOUBLE  *gaussPoints = NULL;
     
     try{
          CALLO_APLCPGN(_mesh->getName(), _newCoordinates->getName(), _zones[i]->getName(),
                       pair_method, &pair_tole, &nbCellMaster, eleMaster.data(), &nbCellSlave, 
                       eleSlave.data(), NodesMaster.data(), &nbNodeMaster, &nb_pairs, &pairs, 
                       &nbInterPoints, &InterSlavePoints, &InterMasterPoints, &gaussPoints);
     }catch(const std::exception& e){
           std::cout << e.what() << std::endl;
     }
     
     // clear 
     this->clear(i);

     // fill the pairing quantities 
     _nbPairs[i] = nb_pairs;
     _listOfPairs[i] = VectorLong(pairs, pairs + 2*nb_pairs );
     _nbIntersectionPoints[i] = VectorLong(nbInterPoints, nbInterPoints + nb_pairs );
     _slaveIntersectionPoints[i] = VectorReal(InterSlavePoints, InterSlavePoints + 16*nb_pairs);
     _masterIntersectionPoints[i] = VectorReal(InterMasterPoints, InterMasterPoints + 16*nb_pairs);
     _gaussPoints[i]  =  VectorReal(gaussPoints, gaussPoints + 72*nb_pairs);

     // free temporary quantities
     free(pairs);
     free(nbInterPoints);
     free(InterSlavePoints);
     free(InterMasterPoints);
     free(gaussPoints);

     return true;
}


ASTERBOOL ContactPairing::clear( ASTERINTEGER i ){

     _listOfPairs[i].clear();
     _nbIntersectionPoints[i].clear();
     _slaveIntersectionPoints[i].clear();
     _masterIntersectionPoints[i].clear();
     _gaussPoints[i].clear();

     return true;
}
