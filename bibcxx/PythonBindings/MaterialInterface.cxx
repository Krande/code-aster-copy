/**
 * @file MaterialInterface.cxx
 * @brief Interface python de Material
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

#include "PythonBindings/MaterialInterface.h"

#include "aster_pybind.h"

void exportMaterialToPython( py::module_ &mod ) {

    py::class_< Material, Material::MaterialPtr, DataStructure >( mod, "Material" )
        .def( py::init( &initFactoryPtr< Material > ) )
        .def( py::init( &initFactoryPtr< Material, std::string > ) )
        .def( py::init( &initFactoryPtr< Material, Material > ) )
        .def( py::init( &initFactoryPtr< Material, Material, VectorString > ) )
        .def( "size", &Material::size, R"(
Return the number of material names.

Returns:
    int: Number of material names.
        )" )
        .def( "getMaterialNames", &Material::getMaterialNames, R"(
Return the list of the material names.

Returns:
    list[str]: List of material names (without "_FO")
        )" )
        .def( "getValueReal", &Material::getValueReal, R"(
Return the value of a property stored as a real.

Raise an exception if the property does not exist.

Arguments:
    materialName (str): Material name (without "_FO").
    propertyName (str): Property name.

Returns:
    float: Property value.
        )",
              py::arg( "materialName" ), py::arg( "propertyName" ) )
        .def( "getValueComplex", &Material::getValueComplex, R"(
Return the value of a property stored as a complex.

Raise an exception if the property does not exist.

Arguments:
    materialName (str): Material name (without "_FO").
    propertyName (str): Property name.

Returns:
    complex: Property value.
        )",
              py::arg( "materialName" ), py::arg( "propertyName" ) )
        .def( "getFunction", &Material::getFunction, R"(
Return the value of a property stored as a function.

Raise an exception if the property does not exist.

Arguments:
    materialName (str): Material name (without "_FO").
    propertyName (str): Property name.

Returns:
    *Function*: Function object, *None* if the property does not exist or is not a function.
        )",
              py::arg( "materialName" ), py::arg( "propertyName" ) )

        .def( "_addProperties", &Material::_addProperties )
        .def( "_storeListReal", &Material::_storeListReal )
        .def( "_storeListFunc", &Material::_storeListFunc )
        .def( "_setTractionFunction", &Material::_setTractionFunction );
};
