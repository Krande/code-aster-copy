/**
 * @file DebugInterface.cxx
 * @brief Python interface for debugging
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

#include "astercxx.h"

#include "PythonBindings/DebugInterface.h"

#include "aster_fort_jeveux.h"
#include "aster_pybind.h"

#include "LinearAlgebra/AssemblyMatrix.h"
#include "LinearAlgebra/ElementaryMatrix.h"
#include "MemoryManager/JeveuxUtils.h"
#include "Meshes/Mesh.h"
#include "Modeling/Model.h"
#include "Numbering/DOFNumbering.h"

#include <string>

static void libaster_debugJeveuxContent( const std::string message ) {
    ASTERINTEGER unit_out = 6;
    std::string base( "G" );
    CALLO_JEIMPR( &unit_out, base, message );
};

void exportDebugToPython( py::module_ &mod ) {

    mod.def( "debugJeveuxContent", &libaster_debugJeveuxContent );
    mod.def( "debugJeveuxExists", &jeveuxExists );
    mod.def( "use_count", &libaster_debugRefCount< MeshPtr > );
    mod.def( "use_count", &libaster_debugRefCount< ModelPtr > );
    mod.def( "use_count", &libaster_debugRefCount< DOFNumberingPtr > );
    mod.def( "use_count", &libaster_debugRefCount< ElementaryMatrixDisplacementRealPtr > );
    mod.def( "use_count", &libaster_debugRefCount< ElementaryMatrixDisplacementComplexPtr > );
    mod.def( "use_count", &libaster_debugRefCount< ElementaryMatrixTemperatureRealPtr > );
    mod.def( "use_count", &libaster_debugRefCount< ElementaryMatrixPressureComplexPtr > );
    mod.def( "use_count", &libaster_debugRefCount< AssemblyMatrixDisplacementRealPtr > );
    mod.def( "use_count", &libaster_debugRefCount< AssemblyMatrixDisplacementComplexPtr > );
    mod.def( "use_count", &libaster_debugRefCount< AssemblyMatrixTemperatureRealPtr > );
    mod.def( "use_count", &libaster_debugRefCount< AssemblyMatrixTemperatureComplexPtr > );
    mod.def( "use_count", &libaster_debugRefCount< AssemblyMatrixPressureRealPtr > );
    mod.def( "use_count", &libaster_debugRefCount< AssemblyMatrixPressureComplexPtr > );
};
