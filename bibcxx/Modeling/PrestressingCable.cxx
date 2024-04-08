/**
 * @file PrestressingCable.cxx
 * @brief Implementation de PrestressingCable
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 22023 EDF R&D                www.code-aster.org
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "astercxx.h"

#include "Modeling/PrestressingCable.h"

#include "Supervis/ResultNaming.h"

PrestressingCable::PrestressingCable( const std::string name, const ModelPtr &model,
                                      const MaterialFieldPtr &mater,
                                      const ElementaryCharacteristicsPtr &cara )
    : DataStructure( name, 8, "CABL_PRECONT" ),
      _model( model ),
      _mater( mater ),
      _cara( cara ),
      _mesh( std::static_pointer_cast< Mesh >( _model->getMesh() ) ),
      _sigin( new ConstantFieldOnCellsReal( getName() + ".CHME.SIGIN", _mesh ) ),
      _cableBP( new Table( getName() + "CABLEBP    " ) ),
      _cableGL( new Table( getName() + "CABLEGL    " ) ),
      _lirela( new ListOfLinearRelationsReal( getName() + ".LIRELA    " ) ),
      _isBuilt( false ) {}

PrestressingCable::PrestressingCable( const ModelPtr &model, const MaterialFieldPtr &mater,
                                      const ElementaryCharacteristicsPtr &cara )
    : PrestressingCable::PrestressingCable( ResultNaming::getNewResultName(), model, mater, cara ) {
}

bool PrestressingCable::build() {
    _cableBP->build();
    _cableGL->build();
    _isBuilt = true;
    return true;
};
