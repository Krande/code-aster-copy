/**
 * @file GeneralizedDOFNumberingInterface.cxx
 * @brief Interface python de GeneralizedDOFNumbering
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "PythonBindings/GeneralizedDOFNumberingInterface.h"

#include "aster_pybind.h"

#include "PythonBindings/VariantModalBasisInterface.h"

void exportGeneralizedDOFNumberingToPython( py::module_ &mod ) {

    py::class_< GeneralizedDOFNumbering, GeneralizedDOFNumbering::GeneralizedDOFNumberingPtr,
                DataStructure >( mod, "GeneralizedDOFNumbering" )
        .def( py::init( &initFactoryPtr< GeneralizedDOFNumbering > ) )
        .def( py::init( &initFactoryPtr< GeneralizedDOFNumbering, std::string > ) )
        .def( "getGeneralizedModel", &GeneralizedDOFNumbering::getGeneralizedModel )
        .def( "getModalBasis", &getModalBasis< GeneralizedDOFNumberingPtr > )
        .def( "setGeneralizedModel", &GeneralizedDOFNumbering::setGeneralizedModel )
        .def( "setModalBasis", py::overload_cast< const ModeResultPtr & >(
                                   &GeneralizedDOFNumbering::setModalBasis ) )
        .def( "setModalBasis", py::overload_cast< const GeneralizedModeResultPtr & >(
                                   &GeneralizedDOFNumbering::setModalBasis ) )
        .def( "getMorseStorage", &GeneralizedDOFNumbering::getMorseStorage );
};
