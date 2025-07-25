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
subroutine testvois(jv_geom, elem_slav_type, &
                    elem_mast_coor, elem_mast_code, elem_slav_nume, &
                    pair_tole, inte_weight, v_mesh_connex, &
                    v_connex_lcum)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/apcoor.h"
#include "asterfort/aptype.h"
#include "asterfort/prjint.h"
#include "asterfort/dctest.h"
!
    integer(kind=8), intent(in) :: jv_geom
    character(len=8), intent(in) :: elem_slav_type
    real(kind=8), intent(in) :: elem_mast_coor(27)
    character(len=8), intent(in) :: elem_mast_code
    integer(kind=8), intent(in) :: elem_slav_nume
    real(kind=8), intent(in) :: pair_tole
    real(kind=8), intent(out) :: inte_weight
    integer(kind=8), pointer :: v_mesh_connex(:)
    integer(kind=8), pointer :: v_connex_lcum(:)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing segment to segment
!
! Compute weight of intersection between slave/master element
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  jv_geom          : JEVEUX adress to updated geometry
! In  elem_slav_type   : type of slave element
! In  elem_mast_coor   : coordinates of master element
! In  elem_mast_code   : code of master element
! In  elem_slav_nume   : index of slave element
! In  pair_tole        : tolerance for pairing
! Out inte_weight      : weight of intersection
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_elin_mast, i_node, i_elin_slav, i_dime
    integer(kind=8) :: elem_slav_nbnode, elem_dime
    real(kind=8) :: elem_slav_coor(27)
    character(len=8) :: elem_slav_code
    integer(kind=8) :: elin_slav_sub(1, 4), elin_mast_sub(1, 4)
    integer(kind=8) :: elin_slav_nbnode(1), elin_mast_nbnode(1)
    integer(kind=8) :: elin_slav_nbsub, elin_mast_nbsub
    real(kind=8) :: elin_slav_coor(27), elin_mast_coor(27)
    character(len=8) :: elin_slav_code, elin_mast_code
    integer(kind=8) :: nb_poin_inte
    real(kind=8) :: ints_weight
    real(kind=8) :: poin_inte(32)
!
! --------------------------------------------------------------------------------------------------
!
    if (elem_slav_nume .ne. 0) then
!
! ----- Get informations about slave element
!
        call aptype(elem_slav_type, &
                    elem_slav_nbnode, elem_slav_code, elem_dime)
!
! ----- Get coordinates of slave element
!
        call apcoor(v_mesh_connex, v_connex_lcum, jv_geom, &
                    elem_slav_nume, elem_slav_nbnode, elem_dime, &
                    elem_slav_coor)
!
! ----- Cut slave element in linearized sub-elements
!
        call dctest(elem_slav_code, elin_slav_sub, elin_slav_nbnode, elin_slav_nbsub, &
                    elin_slav_code)
!
! ----- Cut master element in linearized sub-elements
!
        call dctest(elem_mast_code, elin_mast_sub, elin_mast_nbnode, elin_mast_nbsub, &
                    elin_mast_code)
!
        inte_weight = 0.d0
!
! ----- Loop on linearized slave sub-elements
!
        do i_elin_slav = 1, elin_slav_nbsub
            elin_slav_coor(:) = 0.d0
!
! --------- Get coordinates for current linearized slave sub-element
!
            do i_node = 1, elin_slav_nbnode(i_elin_slav)
                do i_dime = 1, elem_dime
                    elin_slav_coor(3*(i_node-1)+i_dime) = &
                        elem_slav_coor(3*(elin_slav_sub(i_elin_slav, i_node)-1)+i_dime)
                end do
            end do
!
! --------- Loop on linearized master sub-elements
!
            do i_elin_mast = 1, elin_mast_nbsub
                elin_mast_coor(:) = 0.d0
!
! ------------- Get coordinates for current linearized master sub-element
!
                do i_node = 1, elin_mast_nbnode(i_elin_mast)
                    do i_dime = 1, elem_dime
                        elin_mast_coor(3*(i_node-1)+i_dime) = &
                            elem_mast_coor(3*(elin_mast_sub(i_elin_mast, i_node)-1)+i_dime)
                    end do
                end do
!
! ------------- Projection/intersection of elements in slave parametric space
!
                call prjint(pair_tole, elem_dime, &
                            elin_mast_nbnode(i_elin_mast), elin_mast_coor, elin_mast_code, &
                            elin_slav_nbnode(i_elin_slav), elin_slav_coor, elin_slav_code, &
                            poin_inte, ints_weight, nb_poin_inte)
!
                inte_weight = inte_weight+ints_weight
            end do
        end do
    end if
end subroutine
