/**
 * @file Model.cxx
 * @brief Implementation de Model
 * @author Nicolas Sellenet
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

#include "astercxx.h"

#include "Modeling/Partition.h"

const char *const ModelSplitingMethodNames[nbModelSplitingMethod] = { "CENTRALISE", "SOUS_DOMAINE",
                                                                      "GROUP_ELEM" };
const char *const GraphPartitionerNames[nbGraphPartitioner] = { "SCOTCH", "METIS" };

Partition::Partition( const std::string name )
    : DataStructure( name, 8, "PARTITION" ),
      _prti( getName() + ".PRTI" ),
      _prtk( getName() + ".PRTK" ),
      _nupr( getName() + ".PRTI" ),
      _fdim( getName() + ".FDIM" ),
      _feta( getName() + ".FETA" ),
      _fref( getName() + ".FREF" ) {};

const std::string Partition::getMethod() const {
    if ( !_prtk.exists() )
        return ModelSplitingMethodNames[(int)Centralized];
    else {
        _prtk->updateValuePointer();
        return strip( ( *_prtk )[0].toString() );
    }
}
