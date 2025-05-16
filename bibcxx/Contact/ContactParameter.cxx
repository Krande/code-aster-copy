/**
 * @file ContactParameter.cxx
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

#include "Contact/ContactParameter.h"

ContactParameter::ContactParameter( const py::tuple &tup ) : ContactParameter() {
    _algo = tup[0].cast< ContactAlgo >();
    _type = tup[1].cast< ContactType >();
    _vari = tup[2].cast< ContactVariant >();
    _coeff = tup[3].cast< ASTERDOUBLE >();
    _jacType = tup[4].cast< JacobianType >();
};
py::tuple ContactParameter::_getState() const {
    return py::make_tuple( _algo, _type, _vari, _coeff, _jacType );
};

FrictionParameter::FrictionParameter( const py::tuple &tup ) : FrictionParameter() {
    _algo = tup[0].cast< FrictionAlgo >();
    _type = tup[1].cast< FrictionType >();
    _coeff = tup[2].cast< ASTERDOUBLE >();
    _tresca = tup[3].cast< ASTERDOUBLE >();
    _coulomb = tup[4].cast< ASTERDOUBLE >();
};
py::tuple FrictionParameter::_getState() const {
    return py::make_tuple( _algo, _type, _coeff, _tresca, _coulomb );
};

PairingParameter::PairingParameter( const py::tuple &tup ) : PairingParameter() {
    _algo = tup[0].cast< PairingAlgo >();
    _cont_init = tup[1].cast< InitialState >();
    _dist_ratio = tup[2].cast< ASTERDOUBLE >();
    // _beam = (bool)tup[3].cast< int >();
    // _shell = (bool)tup[4].cast< int >();
    _beam = tup[3].cast< bool >();
    _shell = tup[4].cast< bool >();
    _dist_supp = tup[5].cast< GenericFunctionPtr >();
    _cara = tup[6].cast< ElementaryCharacteristicsPtr >();
};
py::tuple PairingParameter::_getState() const {
    return py::make_tuple( _algo, _cont_init, _dist_ratio, (int)_beam, (int)_shell, _dist_supp,
                           _cara );
};
