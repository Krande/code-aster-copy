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
    int i = -1;
    _algo = tup[++i].cast< ContactAlgo >();
    _type = tup[++i].cast< ContactType >();
    _vari = tup[++i].cast< ContactVariant >();
    _coeff = tup[++i].cast< ASTERDOUBLE >();
    _jacType = tup[++i].cast< JacobianType >();
};
py::tuple ContactParameter::_getState() const {
    return py::make_tuple( _algo, _type, _vari, _coeff, _jacType );
};

FrictionParameter::FrictionParameter( const py::tuple &tup ) : FrictionParameter() {
    int i = -1;
    _friction = tup[++i].cast< bool >();
    _algo = tup[++i].cast< FrictionAlgo >();
    _type = tup[++i].cast< FrictionType >();
    _coeff = tup[++i].cast< ASTERDOUBLE >();
    _tresca = tup[++i].cast< ASTERDOUBLE >();
    _coulomb = tup[++i].cast< ASTERDOUBLE >();
};
py::tuple FrictionParameter::_getState() const {
    return py::make_tuple( _friction, _algo, _type, _coeff, _tresca, _coulomb );
};

PairingParameter::PairingParameter( const py::tuple &tup ) : PairingParameter() {
    int i = -1;
    _algo = tup[++i].cast< PairingAlgo >();
    _cont_init = tup[++i].cast< InitialState >();
    _dist_ratio = tup[++i].cast< ASTERDOUBLE >();
    _beam = tup[++i].cast< bool >();
    _shell = tup[++i].cast< bool >();
    _dist_supp = tup[++i].cast< GenericFunctionPtr >();
    _cara = tup[++i].cast< ElementaryCharacteristicsPtr >();
};
py::tuple PairingParameter::_getState() const {
    return py::make_tuple( _algo, _cont_init, _dist_ratio, (int)_beam, (int)_shell, _dist_supp,
                           _cara );
};
