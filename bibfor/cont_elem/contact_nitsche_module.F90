! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
module contact_nitsche_module
!
use contact_type
!
implicit none
!
private
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/matinv.h"
#include "blas/dgemm.h"
#include "blas/dgemv.h"
#include "blas/dger.h"
#include "contact_module.h"
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Nitsche methods
!
! --------------------------------------------------------------------------------------------------
!
    public :: dofsMapping, remappingVect, remappingMatr
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
  function dofsMapping(geom)
!
    implicit none
!
        integer :: dofsMapping(54)
        type(ContactGeom), intent(in) :: geom
!
! --------------------------------------------------------------------------------------------------
!
!   dofs mapping between face to volume
!
! --------------------------------------------------------------------------------------------------
!
        integer :: i_node, i_elem, nb_dofs_slav, nb_dofs_volu
!
        dofsMapping = 0
!
        do i_node = 1, geom%nb_node_slav
            do i_elem = 1, geom%elem_dime
                dofsMapping((i_node-1)*geom%elem_dime+i_elem) = &
                    (geom%mapVolu2Surf(i_node)-1)*geom%elem_dime+i_elem
            end do
        end do
!
        nb_dofs_slav = geom%nb_node_slav * geom%elem_dime
        nb_dofs_volu = geom%nb_node_volu * geom%elem_dime
        do i_node = 1, geom%nb_node_mast
            do i_elem = 1, geom%elem_dime
                dofsMapping(nb_dofs_slav+(i_node-1)*geom%elem_dime+i_elem) = &
                    nb_dofs_volu + (i_node-1)*geom%elem_dime+i_elem
            end do
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine remappingVect(geom, dofsmap, vect, new_vect, coeff)
!
    implicit none
!
        type(ContactGeom), intent(in) :: geom
        integer, intent(in) :: dofsMap(54)
        real(kind=8), intent(in) :: vect(MAX_LAGA_DOFS), coeff
        real(kind=8), intent(inout) :: new_vect(MAX_NITS_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
!   Remapping vector (+= operation)
!
! --------------------------------------------------------------------------------------------------
!
        integer :: i_dof, nb_dofs
!
        nb_dofs = (geom%nb_node_slav + geom%nb_node_mast) * geom%elem_dime
!
        do i_dof = 1, nb_dofs
            new_vect(dofsMap(i_dof)) = new_vect(dofsMap(i_dof)) + coeff * vect(i_dof)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine remappingMatr(geom, dofsmap, matr, new_matr, coeff)
!
    implicit none
!
        type(ContactGeom), intent(in) :: geom
        integer, intent(in) :: dofsMap(54)
        real(kind=8), intent(in) :: matr(MAX_LAGA_DOFS, MAX_LAGA_DOFS), coeff
        real(kind=8), intent(inout) :: new_matr(MAX_NITS_DOFS, MAX_NITS_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
!   Remapping matrix (+= operation)
!
! --------------------------------------------------------------------------------------------------
!
        integer :: i_dof, nb_dofs, j_dof, i_glo, j_glo
!
        nb_dofs = (geom%nb_node_slav + geom%nb_node_mast) * geom%elem_dime
!
        do j_dof = 1, nb_dofs
            j_glo = dofsMap(j_dof)
            do i_dof = 1, nb_dofs
                i_glo = dofsMap(i_dof)
                new_matr(i_glo, j_glo) = new_matr(i_glo, j_glo) + coeff * matr(i_dof, j_dof)
            end do
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
end module
