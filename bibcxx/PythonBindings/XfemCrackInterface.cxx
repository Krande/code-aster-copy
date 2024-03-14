/**
 * @file XfemCrackInterface.cxx
 * @brief Interface python de XfemCrack
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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

#include "PythonBindings/XfemCrackInterface.h"

#include "aster_pybind.h"

void exportXfemCrackToPython( py::module_ &mod ) {

    py::class_< XfemCrack, XfemCrack::XfemCrackPtr, DataStructure >( mod, "XfemCrack" )
        .def( py::init( &initFactoryPtr< XfemCrack, MeshPtr > ) )
        .def( py::init( &initFactoryPtr< XfemCrack, std::string, MeshPtr > ) )
        .def( "build", &XfemCrack::build )
        .def( "updateInternalState", &XfemCrack::update_tables )
        .def( "enrichModelWithXfem", &XfemCrack::enrichModelWithXfem )
        .def( "getMesh", &XfemCrack::getMesh )
        .def( "setMesh", &XfemCrack::setMesh )
        .def( "getAuxiliaryGrid", &XfemCrack::getAuxiliaryGrid )
        .def( "setAuxiliaryGrid", &XfemCrack::setAuxiliaryGrid )
        .def( "getExistingCrackWithGrid", &XfemCrack::getExistingCrackWithGrid )
        .def( "setExistingCrackWithGrid", &XfemCrack::setExistingCrackWithGrid )
        .def( "getDiscontinuityType", &XfemCrack::getDiscontinuityType )
        .def( "setDiscontinuityType", &XfemCrack::setDiscontinuityType )
        .def( "getCrackLipsEntity", &XfemCrack::getCrackLipsEntity )
        .def( "setCrackLipsEntity", &XfemCrack::setCrackLipsEntity )
        .def( "getCrackTipEntity", &XfemCrack::getCrackTipEntity )
        .def( "setCrackTipEntity", &XfemCrack::setCrackTipEntity )
        .def( "getCohesiveCrackTipForPropagation", &XfemCrack::getCohesiveCrackTipForPropagation )
        .def( "setCohesiveCrackTipForPropagation", &XfemCrack::setCohesiveCrackTipForPropagation )
        .def( "getNormalLevelSetFunction", &XfemCrack::getNormalLevelSetFunction )
        .def( "setNormalLevelSetFunction", &XfemCrack::setNormalLevelSetFunction )
        .def( "getTangentialLevelSetFunction", &XfemCrack::getTangentialLevelSetFunction )
        .def( "setTangentialLevelSetFunction", &XfemCrack::setTangentialLevelSetFunction )
        .def( "getCrackShape", &XfemCrack::getCrackShape )
        .def( "setCrackShape", &XfemCrack::setCrackShape )
        .def( "getNormalLevelSetField", &XfemCrack::getNormalLevelSetField )
        .def( "setNormalLevelSetField", &XfemCrack::setNormalLevelSetField )
        .def( "getTangentialLevelSetField", &XfemCrack::getTangentialLevelSetField )
        .def( "setTangentialLevelSetField", &XfemCrack::setTangentialLevelSetField )
        .def( "getEnrichedCells", &XfemCrack::getEnrichedCells )
        .def( "setEnrichedCells", &XfemCrack::setEnrichedCells )
        .def( "getDiscontinuousField", &XfemCrack::getDiscontinuousField )
        .def( "setDiscontinuousField", &XfemCrack::setDiscontinuousField )
        .def( "getEnrichmentType", &XfemCrack::getEnrichmentType )
        .def( "setEnrichmentType", &XfemCrack::setEnrichmentType )
        .def( "getEnrichmentRadiusZone", &XfemCrack::getEnrichmentRadiusZone )
        .def( "setEnrichmentRadiusZone", &XfemCrack::setEnrichmentRadiusZone )
        .def( "getEnrichedLayersNumber", &XfemCrack::getEnrichedLayersNumber )
        .def( "setEnrichedLayersNumber", &XfemCrack::setEnrichedLayersNumber )
        .def( "getJunctingCracks", &XfemCrack::getJunctingCracks )
        .def( "insertJunctingCracks", &XfemCrack::insertJunctingCracks )
        .def( "setPointForJunction", &XfemCrack::setPointForJunction )
        .def( "getCrackTipCoords", &XfemCrack::getCrackTipCoords )
        .def( "getCrackTipBasis", &XfemCrack::getCrackTipBasis )
        .def( "getCrackTipMultiplicity", &XfemCrack::getCrackTipMultiplicity )
        .def( "getTipType", &XfemCrack::getTipType )
        .def( "getCrackTipNodeFacesField", &XfemCrack::getCrackTipNodeFacesField )
        .def( "getCrackFrontRadius", &XfemCrack::getCrackFrontRadius )

        .def( "getTable", &ListOfTables::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              py::arg( "identifier" ) );
};
