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
subroutine lctppe(side, l_axis, l_upda_jaco, &
                  nb_node, elem_dime, elem_code, &
                  elem_init, elem_coor, &
                  gauss_coor, shape_func, &
                  jacobian, norm_g)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/mmnonf.h"
#include "asterfort/mmdonf.h"
#include "asterfort/mmmjac.h"
#include "asterfort/mmnorm.h"
#include "asterfort/mmtang.h"
!
    character(len=*), intent(in) :: side
    integer(kind=8), intent(in) :: elem_dime
    aster_logical, intent(in) :: l_axis, l_upda_jaco
    integer(kind=8), intent(in) :: nb_node
    real(kind=8), intent(in) :: elem_init(nb_node, elem_dime)
    real(kind=8), intent(in) :: elem_coor(nb_node, elem_dime)
    character(len=8), intent(in) :: elem_code
    real(kind=8), intent(in) :: gauss_coor(2)
    real(kind=8), intent(out) :: shape_func(9)
    real(kind=8), intent(out) :: jacobian
    real(kind=8), intent(out) :: norm_g(3)
!
! --------------------------------------------------------------------------------------------------
!
! Contact (LAC) - Elementary computations
!
! Compute geometric quantities for contact matrix
!
! --------------------------------------------------------------------------------------------------
!
! In  side             : side (master or slave)
! In  elem_dime        : dimension of elements
! In  l_axis           : .true. for axisymmetric element
! In  l_upda_jaco      : .true. to use updated jacobian
! In  nb_node          : number of nodes
! In  elem_init        : initial coordinates of element
! In  elem_coor        : updated coordinates of element
! In  elem_code        : code of element
! In  gauss_coor       : coordinates of current integration point
! Out shape_func       : shape functions at integration point
! Out jacobian         : jacobian at integration point
! Out norm_g           : normal at integration point
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_dime, i_node
    real(kind=8) :: elem_coot(3, 9)
    real(kind=8) :: tau1(3), tau2(3), shape_dfunc(2, 9)
!
! --------------------------------------------------------------------------------------------------
!
    shape_func(:) = 0.d0
    shape_dfunc(:, :) = 0.d0
    norm_g(:) = 0.d0
    tau1(:) = 0.d0
    tau2(:) = 0.d0
    jacobian = 0.d0
!
! - Change shape of vector for coordinates of element
!
    elem_coot = 0.d0
    do i_dime = 1, elem_dime
        do i_node = 1, nb_node
            elem_coot(i_dime, i_node) = elem_coor(i_node, i_dime)
        end do
    end do
!
! - Get shape functions and first derivative only (for perf)
!
    call mmnonf(elem_dime, nb_node, elem_code, &
                gauss_coor(1), gauss_coor(2), &
                shape_func)
    call mmdonf(elem_dime, nb_node, elem_code, &
                gauss_coor(1), gauss_coor(2), &
                shape_dfunc)
!
! - Compute normal at integration point
!
    call mmtang(elem_dime, nb_node, elem_coot, shape_dfunc, tau1, tau2)
    call mmnorm(elem_dime, tau1, tau2, norm_g)
    if (side .eq. 'Master') then
        norm_g(1:3) = -norm_g(1:3)
    end if
!
! - Compute updated surfacic jacobian of element
!
    if (l_upda_jaco) then
        call mmmjac(l_axis, nb_node, elem_dime, &
                    elem_code, elem_coot, &
                    shape_func, shape_dfunc, &
                    jacobian)
    else
        elem_coot = 0.d0
        do i_dime = 1, elem_dime
            do i_node = 1, nb_node
                elem_coot(i_dime, i_node) = elem_init(i_node, i_dime)
            end do
        end do
        call mmmjac(l_axis, nb_node, elem_dime, &
                    elem_code, elem_coot, &
                    shape_func, shape_dfunc, &
                    jacobian)
    end if
!
end subroutine
