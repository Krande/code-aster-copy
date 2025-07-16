! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
! ==================================================================================================
!
! Module for types of mesh
!
! ==================================================================================================
module mesh_type
! ==================================================================================================
! ==================================================================================================
    implicit none
! ==================================================================================================
    private
#include "asterf_types.h"
#include "MeshTypes_type.h"
! ==================================================================================================
! Global variables
! ==================================================================================================
! - None
! ==================================================================================================
! Define types
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
! Geometry of a cell
! --------------------------------------------------------------------------------------------------
    type CELL_GEOM
! ----- String for type of cell
        character(len=8) :: cellCode = " "
        aster_logical :: isLinear = ASTER_TRUE
        aster_logical :: isSkin = ASTER_FALSE
! ----- Topological dimension of cell (-1, 0, 1 or 2)
        integer(kind=8) :: cellDime = 0
! ----- Number of nodes
        integer(kind=8) :: nbNode = 0
! ----- Number of neighbours
        integer(kind=8) :: nbNeigh = 0
! ----- Diameter of cell
        real(kind=8) :: diameter = -1.d0
! ----- Barycenter of cell in global basis (X,Y,Z)
        real(kind=8) :: baryGlob(3) = 0.d0
! ----- Barycenter of cell in parametric basis (ksi, eta, zeta)
        real(kind=8) :: baryPara(3) = 0.d0
! ----- Coordinates of cell in global basis (X,Y,Z)
        real(kind=8), dimension(3, MT_NNOMAX3D) :: coorNodeGlob = 0.d0
! ----- Coordinates of cell in parametric basis (ksi, eta, zeta)
        real(kind=8), dimension(3, MT_NNOMAX3D) :: coorNodePara = 0.d0
    end type CELL_GEOM
! --------------------------------------------------------------------------------------------------
! Base for skin cell
! --------------------------------------------------------------------------------------------------
    type CELL_SKIN_BASE
! ----- Dimension of space
        integer(kind=8) :: spaceDime = 0
! ----- Flag for _external_ normal
        aster_logical :: normIsExte = ASTER_TRUE
! ----- Normal
        real(kind=8) :: norm(3) = 0.d0
! ----- Tangents
        real(kind=8) :: tau(3, 2) = 0.d0
    end type CELL_SKIN_BASE
!===================================================================================================
!===================================================================================================
    public :: CELL_GEOM, CELL_SKIN_BASE
contains
!===================================================================================================
!===================================================================================================
end module mesh_type
