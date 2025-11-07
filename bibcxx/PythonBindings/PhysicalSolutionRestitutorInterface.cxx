/**
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

#include "PythonBindings/PhysicalSolutionRestitutorInterface.h"

#include "aster_pybind.h"

void exportPhysicalSolutionRestitutorToPython( py::module_ &mod ) {
    py::class_< PhysicalSolutionRestitutor, std::shared_ptr< PhysicalSolutionRestitutor > >(
        mod, "PhysicalSolutionRestitutor" )

        .def( py::init( []( const TransientGeneralizedResultPtr &resgen ) {
                  return initFactoryPtr< PhysicalSolutionRestitutor, TransientGeneralizedResultPtr,
                                         ASTERINTEGER >( resgen, 1, 0 );
              } ),
              R"(
PhysicalSolutionRestitutor(resgen)

Create a PhysicalSolutionRestitutor associated with a transient generalized result.

Arguments:
    resgen (TransientGeneralizedResult):
        Pointer to the transient generalized result from which physical fields
        will be reconstructed.

Returns:
    PhysicalSolutionRestitutor:
        A restitution object initialized with the given transient result.
              )",
              py::arg( "resgen" ) )

        .def( py::init( &initFactoryPtr< PhysicalSolutionRestitutor, TransientGeneralizedResultPtr,
                                         ASTERINTEGER, ASTERINTEGER > ),
              py::arg( "resgen" ), py::arg( "nbatch" ), py::arg( "ar" ),
              R"(
PhysicalSolutionRestitutor(resgen, nbatch, ar)

Create a PhysicalSolutionRestitutor with control over batching and activation parameters.

Arguments:
    resgen (TransientGeneralizedResult):
        Transient generalized result object containing modal information.
    nbatch (int):
        Number of batches used during restitution computations.
    ar (int):
        Additional parameter for enveloppe computations.

Returns:
    PhysicalSolutionRestitutor:
        Configured restitution object for field reconstruction.
              )" )

        .def( "get_displacement_coeffs", &PhysicalSolutionRestitutor::getDisplacementCoeffs,
              py::return_value_policy::reference_internal,
              R"(
get_displacement_coeffs() -> list[float]

Return the generalized displacement coefficients.

These coefficients represent the modal displacement amplitudes
used during physical restitution.

Returns:
    list[float]:
        Reference to the vector of displacement coefficients (no copy).
              )" )

        .def( "get_velocity_coeffs", &PhysicalSolutionRestitutor::getVelocityCoeffs,
              py::return_value_policy::reference_internal,
              R"(
get_velocity_coeffs() -> list[float]

Return the generalized velocity coefficients.

These coefficients represent the modal velocity amplitudes
used to reconstruct physical velocity fields.

Returns:
    list[float]:
        Reference to the vector of velocity coefficients (no copy).
              )" )

        .def( "get_acceleration_coeffs", &PhysicalSolutionRestitutor::getAccelerationCoeffs,
              py::return_value_policy::reference_internal,
              R"(
get_acceleration_coeffs() -> list[float]

Return the generalized acceleration coefficients.

These coefficients represent the modal acceleration amplitudes
used to reconstruct physical acceleration fields.

Returns:
    list[float]:
        Reference to the vector of acceleration coefficients (no copy).
              )" )

        .def( "computeMaxForFieldsOnNodes", &PhysicalSolutionRestitutor::computeMaxForFieldsOnNodes,
              R"(
computeMaxForFieldsOnNodes() -> dict[str, FieldOnNodesReal]

Compute the time-maximum of all modal fields defined on nodes.

This function processes all nodal fields associated with the transient result
and returns the maximum (component-wise or field-wise) observed over time
for each field.

Returns:
    dict[str, FieldOnNodesReal]:
        Mapping between field names and their corresponding 
        nodal fields containing the maximum values over time.
              )" )

        .def( "computeMaxForFieldsOnCells", &PhysicalSolutionRestitutor::computeMaxForFieldsOnCells,
              R"(
computeMaxForFieldsOnCells() -> dict[str, FieldOnCellsReal]

Compute the time-maximum of all modal fields defined on cells.

Similar to `computeMaxForFieldsOnNodes`, but applied to cell-based fields.
Each entry in the returned dictionary corresponds to a field name and its
cell field containing maximum values over the transient duration.

Returns:
    dict[str, FieldOnCellsReal]:
        Mapping between field names and their corresponding 
        cell fields containing the maximum values over time.
              )" );
}
