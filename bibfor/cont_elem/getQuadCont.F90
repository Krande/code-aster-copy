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
subroutine getQuadCont(elem_dime, l_axis, nb_node_slav, &
                        elem_slav_code, elem_slav_coor,&
                        elem_mast_code, &
                        nb_qp, coor_qp, weight_qp )
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/lcptga.h"
#include "asterfort/latrco.h"
#include "asterfort/jevech.h"
#include "asterfort/mmmjac.h"
#include "asterfort/mmnonf.h"
#include "asterfort/mmdonf.h"
#include "jeveux.h"
!
integer, intent(in) :: elem_dime, nb_node_slav
aster_logical, intent(in) :: l_axis
character(len=8), intent(in) :: elem_slav_code, elem_mast_code
real(kind=8), intent(in) :: elem_slav_coor(3, 9)
real(kind=8), intent(out) :: coor_qp(2, 48), weight_qp(48)
integer, intent(out) :: nb_qp
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Quadrature
!
! Compute quadrature in slave slide
!
! --------------------------------------------------------------------------------------------------
!
! In  elem_dime        : dimension of elements
! In  elem_slav_code   : code element for slave side from contact element
! In  elem_slav_coor   : coordinates from slave side of contact element
! Out  nb_qp           : number of quadrature points
! Out  coor_qp         : coordinates of quadrature points (parametric slave space)
! Out  weight_qp       : weight of quadrature points
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_slav_line, l_mast_line
    integer :: i_node, i_dime, i_tria, i_gauss, jcont
    integer :: nb_tria, nb_poin_inte, nb_gauss
    real(kind=8) :: tria_coot_sl(2,3), tria_coor_sl(16), poin_inte_sl(16)
    real(kind=8) :: gauss_weight_sl(12), gauss_coor_sl(2,12)
    real(kind=8) :: shape_func(9), shape_dfunc(2, 9), jacobian_sl
    character(len=8) :: elga_fami
!
! --------------------------------------------------------------------------------------------------
!
    nb_qp = 0
    coor_qp = 0
    weight_qp = 0
!
! - Get intersection point
!
    call jevech('PCONFR', 'L', jcont)
    nb_poin_inte = int(zr(jcont-1+1))
    poin_inte_sl = zr(jcont-1+2:jcont-1+17)
!
    l_slav_line = elem_slav_code == "SE2" .or. elem_slav_code == "TR3" .or. elem_slav_code == "QU4"
    l_mast_line = elem_mast_code == "SE2" .or. elem_mast_code == "TR3" .or. elem_mast_code == "QU4"
!
! - Triangulation of convex polygon defined by intersection points
    if (elem_dime .eq. 3) then
        if(nb_poin_inte == 3) then
            nb_tria = 1
        else
            nb_tria = nb_poin_inte
        end if
!
        if(l_slav_line .and. l_mast_line) then
            ! order 3 by triangle
            elga_fami = 'FPG4'
        else
            ! order 5 by triangle
            elga_fami = 'FPG7'
        end if
    elseif (elem_dime .eq. 2) then
        nb_tria = 1
        if(l_slav_line .and. l_mast_line) then
            ! order 5
            elga_fami = 'FPG3'
        else
            ! order 7
            elga_fami = 'FPG4'
        end if
    else
        ASSERT(ASTER_FALSE)
    end if
! - Loop on triangles
    do i_tria = 1, nb_tria
! ----- Coordinates of current triangle (slave)
        tria_coor_sl(:) = 0.d0
        if (elem_dime .eq. 3) then
            call latrco(i_tria, nb_poin_inte, poin_inte_sl, tria_coor_sl)
        elseif (elem_dime .eq. 2) then
            tria_coor_sl(1:16) = poin_inte_sl(1:16)
        endif
! ----- Change shape of vector (slave)
        tria_coot_sl(1:2,1:3) = 0.d0
        if (elem_dime .eq. 3) then
            do i_node = 1,3
                do i_dime = 1,(elem_dime-1)
                    tria_coot_sl(i_dime, i_node) = &
                        tria_coor_sl((i_node-1)*(elem_dime-1)+i_dime)
                end do
            end do
        else
            tria_coot_sl(1,1) = tria_coor_sl(1)
            tria_coot_sl(2,1) = 0.d0
            tria_coot_sl(1,2) = tria_coor_sl(2)
            tria_coot_sl(2,2) = 0.d0
        end if
! ----- Get integration points for slave element
        call lcptga(elem_dime, tria_coot_sl , elga_fami      ,&
                    nb_gauss , gauss_coor_sl, gauss_weight_sl)
! ----- Loop on integration points
        do i_gauss = 1, nb_gauss
            nb_qp = nb_qp + 1
            ASSERT(nb_qp <= 48)
! --------- Get current integration point (slave)
            do i_dime = 1, elem_dime-1
                coor_qp(i_dime, nb_qp) = gauss_coor_sl(i_dime, i_gauss)
            end do
!
! - Get shape functions and first derivative only (for perf)
!
            call mmnonf(elem_dime, nb_node_slav, elem_slav_code,&
                        gauss_coor_sl(1, i_gauss), gauss_coor_sl(2, i_gauss),&
                        shape_func )
            call mmdonf(elem_dime, nb_node_slav, elem_slav_code,&
                        gauss_coor_sl(1, i_gauss), gauss_coor_sl(2, i_gauss),&
                        shape_dfunc)
! --------- Compute jacobian
            call mmmjac(l_axis, nb_node_slav, elem_dime,&
                        elem_slav_code , elem_slav_coor  ,&
                        shape_func, shape_dfunc, jacobian_sl)

            weight_qp(nb_qp) = jacobian_sl * gauss_weight_sl(i_gauss)
        end do
    end do
!
end subroutine
