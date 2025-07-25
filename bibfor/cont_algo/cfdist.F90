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

subroutine cfdist(ds_contact, i_zone, elem_slav_indx, poin_coor, time_curr, &
                  gap_user, node_slav_indx_)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cfdism.h"
#include "asterfort/fointe.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mminfl.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    type(NL_DS_Contact), intent(in) :: ds_contact
    integer(kind=8), intent(in) :: i_zone
    integer(kind=8), intent(in) :: elem_slav_indx
    real(kind=8), intent(in) :: poin_coor(3)
    real(kind=8), intent(in) :: time_curr
    real(kind=8), intent(out) :: gap_user
    integer(kind=8), optional, intent(in) :: node_slav_indx_
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Continue/Discrete method - Compute user gap
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_contact       : datastructure for contact management
! In  i_zone           : index of contact zone
! In  elem_slav_indx   : index of slave element (in contact datastructure)
! In  time_curr        : current time
! In  poin_coor        : coordinates of (contact) point
! In  node_slav_indx   : index of slave node (in contact datastructure)
! Out gap_user         : user gap (from DIST_* keywords)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ier
    character(len=8) :: para_name(4)
    real(kind=8) :: para_vale(4)
    real(kind=8) :: gap_user_mast, gap_user_slav, gap_structural
    character(len=8) :: gap_mast_func, gap_slav_func
    aster_logical :: l_dist_beam, l_dist_shell, l_dist_slav, l_dist_mast
    character(len=24) :: sdcont_jeucoq
    real(kind=8), pointer :: v_sdcont_jeucoq(:) => null()
    character(len=24) :: sdcont_jeupou
    real(kind=8), pointer :: v_sdcont_jeupou(:) => null()
    character(len=24) :: sdcont_jeufo1
    character(len=8), pointer :: v_sdcont_jeufo1(:) => null()
    character(len=24) :: sdcont_jeufo2
    character(len=8), pointer :: v_sdcont_jeufo2(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    gap_user = 0.d0
    gap_user_mast = 0.d0
    gap_user_slav = 0.d0
    gap_structural = 0.d0
!
! - Acces to contact objects
!
    sdcont_jeucoq = ds_contact%sdcont_defi(1:16)//'.JEUCOQ'
    sdcont_jeupou = ds_contact%sdcont_defi(1:16)//'.JEUPOU'
    sdcont_jeufo1 = ds_contact%sdcont_defi(1:16)//'.JFO1CO'
    sdcont_jeufo2 = ds_contact%sdcont_defi(1:16)//'.JFO2CO'
    call jeveuo(sdcont_jeucoq, 'L', vr=v_sdcont_jeucoq)
    call jeveuo(sdcont_jeupou, 'L', vr=v_sdcont_jeupou)
    call jeveuo(sdcont_jeufo1, 'L', vk8=v_sdcont_jeufo1)
    call jeveuo(sdcont_jeufo2, 'L', vk8=v_sdcont_jeufo2)
!
! - Set parameters for evaluate functions
!
    para_name(1) = 'X'
    para_name(2) = 'Y'
    para_name(3) = 'Z'
    para_name(4) = 'INST'
    para_vale(1) = poin_coor(1)
    para_vale(2) = poin_coor(2)
    para_vale(3) = poin_coor(3)
    para_vale(4) = time_curr
!
! - Supplementary gaps
!
    l_dist_beam = mminfl(ds_contact%sdcont_defi, 'DIST_POUTRE', i_zone)
    l_dist_shell = mminfl(ds_contact%sdcont_defi, 'DIST_COQUE', i_zone)
    l_dist_mast = mminfl(ds_contact%sdcont_defi, 'DIST_MAIT', i_zone)
    l_dist_slav = mminfl(ds_contact%sdcont_defi, 'DIST_ESCL', i_zone)
!
! - Evaluate DIST_MAIT
!
    if (l_dist_mast) then
        gap_mast_func = v_sdcont_jeufo1(i_zone)
        call fointe('F', gap_mast_func, 4, para_name, para_vale, &
                    gap_user_mast, ier)
    end if
!
! - Evaluate DIST_ESCL
!
    if (l_dist_slav) then
        gap_slav_func = v_sdcont_jeufo2(i_zone)
        call fointe('F', gap_slav_func, 4, para_name, para_vale, &
                    gap_user_slav, ier)
    end if
!
! - Evaluate DIST_POUTRE/DIST_COQUE
!
    if (l_dist_shell .or. l_dist_beam) then
        if (present(node_slav_indx_)) then
            call cfdism(ds_contact, l_dist_beam, l_dist_shell, node_slav_indx_, &
                        gap_structural)
        else
            if (l_dist_beam) then
                gap_structural = gap_structural+v_sdcont_jeupou(elem_slav_indx)
            end if
            if (l_dist_shell) then
                gap_structural = gap_structural+v_sdcont_jeucoq(elem_slav_indx)
            end if
        end if
    end if
!
! - Total user gap
!
    gap_user = gap_user_mast+gap_user_slav+gap_structural
!
end subroutine
