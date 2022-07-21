/**
 * @file FiniteElementDescriptorInterface.cxx
 * @brief Interface python de FiniteElementDescriptor
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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

#include "PythonBindings/FiniteElementDescriptorInterface.h"

#include "aster_pybind.h"

void exportFiniteElementDescriptorToPython( py::module_ &mod ) {

    py::class_< FiniteElementDescriptor, FiniteElementDescriptor::FiniteElementDescriptorPtr,
                DataStructure >( mod, "FiniteElementDescriptor" )
        .def( py::init( &initFactoryPtr< FiniteElementDescriptor, BaseMeshPtr > ) )
        .def( py::init( &initFactoryPtr< FiniteElementDescriptor, std::string, BaseMeshPtr > ) )
        .def( py::init(
            &initFactoryPtr< FiniteElementDescriptor, FiniteElementDescriptorPtr, VectorString > ) )
        .def( py::init( &initFactoryPtr< FiniteElementDescriptor, ModelPtr, VectorString > ) )
        .def( "getPhysics", &FiniteElementDescriptor::getPhysics )
        .def( "getMesh", &FiniteElementDescriptor::getMesh )
        .def( "getModel", &FiniteElementDescriptor::getModel )
        .def( "setModel", &FiniteElementDescriptor::setModel )
        .def( "getNumberOfVirtualNodesobj", &FiniteElementDescriptor::getNumberOfVirtualNodesobj )
        .def( "getNema", &FiniteElementDescriptor::getNema )
        .def( "getListOfGroupOfElements", &FiniteElementDescriptor::getListOfGroupOfElements )
#ifdef ASTER_HAVE_MPI
        .def( "transferDofDescriptorFrom", &FiniteElementDescriptor::transferDofDescriptorFrom )
#endif /* ASTER_HAVE_MPI */
        ;
};
