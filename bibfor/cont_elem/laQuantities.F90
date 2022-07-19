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
subroutine laQuantities(elem_dime, nb_node_slav, nb_node_mast, &
                        indi_lagc, &
                        slav_coor_init, mast_coor_init, &
                        slav_coor_curr, mast_coor_curr, &
                        slav_depl_curr, mast_depl_curr, lagc_curr)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "jeveux.h"
#include "Contact_type.h"
!
integer, intent(in) :: elem_dime, indi_lagc(9)
integer, intent(in) :: nb_node_slav, nb_node_mast
real(kind=8), intent(out) :: mast_coor_init(3, 9), slav_coor_init(3, 9)
real(kind=8), intent(out) :: mast_coor_curr(3, 9), slav_coor_curr(3, 9)
real(kind=8), intent(out) :: mast_depl_curr(3, 9), slav_depl_curr(3, 9)
real(kind=8), intent(out) :: lagc_curr(4)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Elementary computations
!
! Get physical quantities
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i_node_slav, i_node_mast, i_dime, nb_lagr
    integer :: jv_geom, jv_disp_incr, jv_disp, jv_geom_c
    real(kind=8) :: mast_depl_incr(3, 9), slav_depl_incr(3, 9)
    real(kind=8) :: mast_depl_prev(3, 9), slav_depl_prev(3, 9)
!
! --------------------------------------------------------------------------------------------------
!
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PGEOMCR', 'L', jv_geom_c)
    call jevech('PDEPL_P', 'L', jv_disp_incr)
    call jevech('PDEPL_M', 'L', jv_disp )
!
! - Initializations
!
    mast_coor_init = 0.d0
    slav_coor_init = 0.d0
    mast_coor_curr = 0.d0
    slav_coor_curr = 0.d0
    mast_depl_prev = 0.d0
    slav_depl_prev = 0.d0
    mast_depl_incr = 0.d0
    slav_depl_incr = 0.d0
    lagc_curr = 0.d0
!
! - Slave nodes
!
    nb_lagr = 0
    do i_node_slav = 1, nb_node_slav
        do i_dime = 1, elem_dime
            slav_coor_init(i_dime, i_node_slav) =&
                zr(jv_geom+(i_node_slav-1)*elem_dime+i_dime-1)
            slav_coor_curr(i_dime, i_node_slav) =&
                zr(jv_geom_c+(i_node_slav-1)*elem_dime+i_dime-1)
            slav_depl_prev(i_dime, i_node_slav) =&
                zr(jv_disp+(i_node_slav-1)*elem_dime+i_dime-1 + nb_lagr)
            slav_depl_incr(i_dime, i_node_slav) =&
                zr(jv_disp_incr+(i_node_slav-1)*elem_dime+i_dime-1 + nb_lagr)
        end do
        if( indi_lagc(i_node_slav) == 1) then
            nb_lagr =  nb_lagr + 1
            lagc_curr(nb_lagr) = zr(jv_disp+(i_node_slav-1)*elem_dime+elem_dime-1 + nb_lagr) + &
                zr(jv_disp_incr+(i_node_slav-1)*elem_dime+elem_dime-1 + nb_lagr)
        end if
    end do
!
! - Master nodes
!
    do i_node_mast = 1, nb_node_mast
        do i_dime = 1, elem_dime
            mast_coor_init(i_dime, i_node_mast) = &
                zr(jv_geom+nb_node_slav*elem_dime+(i_node_mast-1)*elem_dime+i_dime- 1)
            mast_coor_curr(i_dime, i_node_mast) = &
                zr(jv_geom_c+nb_node_slav*elem_dime+(i_node_mast-1)*elem_dime+i_dime- 1)
            mast_depl_prev(i_dime, i_node_mast) = &
                zr(jv_disp+nb_node_slav*elem_dime+nb_lagr+(i_node_mast-1)*elem_dime+i_dime- 1)
            mast_depl_incr(i_dime, i_node_mast) = &
                zr(jv_disp_incr+nb_node_slav*elem_dime+nb_lagr+(i_node_mast-1)*elem_dime+i_dime- 1)
        end do
    end do
!
    slav_depl_curr = slav_depl_prev + slav_depl_incr
    mast_depl_curr = mast_depl_prev + mast_depl_incr
!
end subroutine
