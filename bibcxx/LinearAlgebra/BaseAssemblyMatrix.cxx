/**
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

#include "LinearAlgebra/BaseAssemblyMatrix.h"

#include "aster_fort_calcul.h"

BaseAssemblyMatrix::BaseAssemblyMatrix( const std::string &name, const std::string &type )
    : DataStructure( name, 19, type ),
      _description( JeveuxVectorChar24( getName() + ".REFA" ) ),
      _scaleFactorLagrangian( JeveuxVectorReal( getName() + ".CONL" ) ),
      _listOfElementaryMatrix( JeveuxVectorChar24( getName() + ".LIME" ) ),
      _perm( JeveuxVectorLong( getName() + ".PERM" ) ),
      _ccid( JeveuxVectorLong( getName() + ".CCID" ) ),
      _ccll( JeveuxVectorLong( getName() + ".CCLL" ) ),
      _ccii( JeveuxVectorLong( getName() + ".CCII" ) ),
      _isEmpty( true ),
      _isFactorized( false ),
      _dofNum( nullptr ),
      _listOfLoads( std::make_shared< ListOfLoads >() ){};

BaseAssemblyMatrix::BaseAssemblyMatrix( const PhysicalProblemPtr phys_prob,
                                        const std::string &type )
    : BaseAssemblyMatrix( type ) {
    _dofNum = phys_prob->getDOFNumbering();
    _listOfLoads = phys_prob->getListOfLoads();
};

BaseAssemblyMatrix::BaseAssemblyMatrix( const std::string &name, const std::string &type,
                                        const BaseAssemblyMatrix &toCopy )
    : BaseAssemblyMatrix( name, type ) {
    // Jeveux Pointer
    ( *_description ) = ( *toCopy._description );
    ( *_scaleFactorLagrangian ) = ( *toCopy._scaleFactorLagrangian );
    ( *_listOfElementaryMatrix ) = ( *toCopy._listOfElementaryMatrix );
    ( *_perm ) = ( *toCopy._perm );
    ( *_ccid ) = ( *toCopy._ccid );
    ( *_ccll ) = ( *toCopy._ccll );
    ( *_ccii ) = ( *toCopy._ccii );
    // Objects
    _dofNum = toCopy._dofNum;
    _listOfLoads = toCopy._listOfLoads;
    _isEmpty = toCopy._isEmpty;
    _isFactorized = toCopy._isFactorized;
}

BaseAssemblyMatrix::BaseAssemblyMatrix( BaseAssemblyMatrix &&other )
    : DataStructure{ std::move( other ) } {
    // Jeveux Pointer
    _description = other._description;
    _scaleFactorLagrangian = other._scaleFactorLagrangian;
    _listOfElementaryMatrix = other._listOfElementaryMatrix;
    _perm = other._perm;
    _ccid = other._ccid;
    _ccll = other._ccll;
    _ccii = other._ccii;
    // Objects
    _dofNum = other._dofNum;
    _listOfLoads = other._listOfLoads;
    _isEmpty = other._isEmpty;
    _isFactorized = other._isFactorized;
}

void BaseAssemblyMatrix::updateDOFNumbering(){
    if ( _description->exists() ){
        _description->updateValuePointer();
        std::string dofName = ( *_description )[1];
        _dofNum = std::make_shared< DOFNumbering > ( dofName );
    }
};

void BaseAssemblyMatrix::symmetrize() { CALL_MATR_ASSE_SYME( getName() ); };

bool BaseAssemblyMatrix::isMPIFull() {
    _description->updateValuePointer();
    return trim( ( *_description )[10].toString() ) == "MPI_COMPLET";
};
