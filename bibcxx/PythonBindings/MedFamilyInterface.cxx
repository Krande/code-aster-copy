/**
 * @file MedFamilyInterface.cxx
 * @brief Interface python de MedFamily
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

#include "PythonBindings/MedFamilyInterface.h"

#include "aster_pybind.h"

// Not DataStructures
// aslint: disable=C3006

#ifdef ASTER_HAVE_MED
void exportMedFamilyToPython( py::module_ &mod ) {

    py::class_< MedFamily, MedFamily::MedFamilyPtr >( mod, "MedFamily" )
        .def( "__pickling_disabled__", disable_pickling< MedFamily >() )

        .def( "getGroups", &MedFamily::getGroups, R"(
Get list of groups in family )" )
        .def( "getName", &MedFamily::getName, R"(
Get family name )" )
        .def( "getId", &MedFamily::getId, R"(
Get family med id )" );
};

#endif
