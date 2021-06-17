/**
 * @file ElementaryMatrixInterface.cxx
 * @brief Interface python de ElementaryMatrix
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

#include <boost/python.hpp>

namespace py = boost::python;
#include <PythonBindings/factory.h>
#include "PythonBindings/ElementaryMatrixInterface.h"

void exportElementaryMatrixToPython() {

    py::class_< BaseElementaryMatrix, BaseElementaryMatrixPtr,
            py::bases< DataStructure > >( "ElementaryMatrixDisplacementReal", py::no_init )
    // fake initFactoryPtr: not buildable
    // fake initFactoryPtr: not buildable
        .def( "addFiniteElementDescriptor",
              &ElementaryMatrixDisplacementReal::addFiniteElementDescriptor )
        .def( "getFiniteElementDescriptors",
              &ElementaryMatrixDisplacementReal::getFiniteElementDescriptors )
        .def( "getMaterialField", &ElementaryMatrixDisplacementReal::getMaterialField )
        .def( "getModel", &ElementaryMatrixDisplacementReal::getModel )
        .def( "getMesh", &ElementaryMatrixDisplacementReal::getMesh )
        .def( "setMaterialField", &ElementaryMatrixDisplacementReal::setMaterialField )
        .def( "setModel", &ElementaryMatrixDisplacementReal::setModel );

    py::class_< ElementaryMatrixDisplacementReal,
            ElementaryMatrixDisplacementReal::ElementaryMatrixPtr,
            py::bases< BaseElementaryMatrix > >( "ElementaryMatrixDisplacementReal",
                                                     py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ElementaryMatrixDisplacementReal >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ElementaryMatrixDisplacementReal,
                                                std::string >))
        .def( "build", &ElementaryMatrixDisplacementReal::build );

    py::class_< ElementaryMatrixDisplacementComplex,
            ElementaryMatrixDisplacementComplex::ElementaryMatrixPtr,
            py::bases< BaseElementaryMatrix > >( "ElementaryMatrixDisplacementComplex",
                                                     py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ElementaryMatrixDisplacementComplex >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ElementaryMatrixDisplacementComplex,
                                                std::string >))
        .def( "build", &ElementaryMatrixDisplacementComplex::build );

    py::class_< ElementaryMatrixTemperatureReal,
            ElementaryMatrixTemperatureReal::ElementaryMatrixPtr,
            py::bases< BaseElementaryMatrix > >( "ElementaryMatrixTemperatureReal",
                                                     py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ElementaryMatrixTemperatureReal >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ElementaryMatrixTemperatureReal,
                                                std::string >))
        .def( "build", &ElementaryMatrixTemperatureReal::build );

    py::class_< ElementaryMatrixPressureComplex,
            ElementaryMatrixPressureComplex::ElementaryMatrixPtr,
            py::bases< BaseElementaryMatrix > >( "ElementaryMatrixPressureComplex",
                                                     py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ElementaryMatrixPressureComplex >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ElementaryMatrixPressureComplex,
                                                std::string >))
        .def( "build", &ElementaryMatrixPressureComplex::build );
};
