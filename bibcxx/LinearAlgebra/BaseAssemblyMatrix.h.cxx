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

BaseAssemblyMatrix::BaseAssemblyMatrix( const std::string &name, const std::string &type )
    : DataStructure( name, 19, type ), _description( JeveuxVectorChar24( getName() + ".REFA" ) ),
      _scaleFactorLagrangian( JeveuxVectorReal( getName() + ".CONL" ) ),
      _listOfElementaryMatrix( JeveuxVectorChar24( getName() + ".LIME" ) ),
      _perm( JeveuxVectorLong( getName() + ".PERM" ) ),
      _ccid( JeveuxVectorLong( getName() + ".CCID" ) ),
      _ccll( JeveuxVectorLong( getName() + ".CCLL" ) ),
      _ccii( JeveuxVectorLong( getName() + ".CCII" ) ),
      _isEmpty( true ),
      _isFactorized( false ),
      _dofNum( nullptr ),
      _listOfLoads( boost::make_shared< ListOfLoads >() ){};

BaseAssemblyMatrix::BaseAssemblyMatrix( const PhysicalProblemPtr phys_prob,
                                        const std::string &type )
    : BaseAssemblyMatrix( type ) {
    _dofNum = phys_prob->getDOFNumbering();
    _listOfLoads = phys_prob->getListOfLoads();
};
