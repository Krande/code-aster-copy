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
subroutine laQuantities(geom)
!
use contact_module
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "jeveux.h"
!
type(ContactGeom), intent(inout) :: geom
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Elementary computations
!
! Get physical quantities
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_fric
    integer :: i_node_slav, i_node_mast, i_dime, nb_lagr, nb_lagr_c, elem_dime, nb_node_slav
    integer :: jv_geom, jv_disp_incr, jv_disp, jv_geom_c, index
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
    mast_depl_prev = 0.d0
    slav_depl_prev = 0.d0
    mast_depl_incr = 0.d0
    slav_depl_incr = 0.d0
!
    l_fric = lteatt('FROTTEMENT','OUI')
    elem_dime = geom%elem_dime
    nb_node_slav = geom%nb_node_slav
!
! - Slave nodes
!
    nb_lagr = 0
    nb_lagr_c = 0
    index = 0
    do i_node_slav = 1, nb_node_slav
        do i_dime = 1, elem_dime
            geom%slav_coor_init(i_dime, i_node_slav) = zr(jv_geom-1+index+i_dime)
            geom%slav_coor_curr(i_dime, i_node_slav) = zr(jv_geom_c-1+index+i_dime)
            slav_depl_prev(i_dime, i_node_slav) = zr(jv_disp-1+index+i_dime + nb_lagr)
            slav_depl_incr(i_dime, i_node_slav) = zr(jv_disp_incr-1+index+i_dime + nb_lagr)
        end do
!
        index = index + elem_dime
!
        if( geom%indi_lagc(i_node_slav) > 0) then
            nb_lagr_c = nb_lagr_c + 1
            geom%slav_lagc_curr(nb_lagr_c) = &
                   zr(jv_disp-1+index + nb_lagr + 1) +  zr(jv_disp_incr-1+index + nb_lagr + 1)
            if(l_fric) then
                geom%slav_lagf_curr(1, nb_lagr_c) = &
                    zr(jv_disp-1+index + nb_lagr + 2) +  zr(jv_disp_incr-1+index + nb_lagr + 2)
                if(geom%elem_dime == 3) then
                    geom%slav_lagf_curr(2, nb_lagr_c) = &
                        zr(jv_disp-1+index + nb_lagr + 3) +  zr(jv_disp_incr-1+index + nb_lagr + 3)
                end if
            end if
            nb_lagr =  nb_lagr + geom%indi_lagc(i_node_slav)
        end if
    end do
!
! - Master nodes
!
    do i_node_mast = 1, geom%nb_node_mast
        do i_dime = 1, elem_dime
            geom%mast_coor_init(i_dime, i_node_mast) = zr(jv_geom-1+index+i_dime)
            geom%mast_coor_curr(i_dime, i_node_mast) = zr(jv_geom_c-1+index+i_dime)
            mast_depl_prev(i_dime, i_node_mast) = zr(jv_disp-1+index+nb_lagr+i_dime)
            mast_depl_incr(i_dime, i_node_mast) = zr(jv_disp_incr-1+index+nb_lagr+i_dime)
        end do
        index = index + elem_dime
    end do
!
    geom%slav_depl_curr = slav_depl_prev + slav_depl_incr
    geom%mast_depl_curr = mast_depl_prev + mast_depl_incr
!
end subroutine
