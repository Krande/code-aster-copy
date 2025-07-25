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

subroutine get_patch_info(sdappa, patch_indx, nb_elem_patch, list_elem)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/jeveuo.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=19), intent(in) :: sdappa
    integer(kind=8), intent(in) :: patch_indx
    integer(kind=8), intent(out) :: nb_elem_patch
    integer(kind=8), intent(out) :: list_elem(5)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing
!
! Segment to segment - Information from patch
!
! --------------------------------------------------------------------------------------------------
!
! In  sdappa           : name of pairing datastructure
! In  patch_indx       : index of patch
! Out nb_elem_patch    : number of elements in patch
! Out list_elem        : list of elements in patch
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_elem_patch
    character(len=24) :: sdappa_info
    integer(kind=8), pointer :: v_sdappa_info(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nb_elem_patch = 0
    list_elem(1:5) = 0
!
! - Access to datastructure
!
    sdappa_info = sdappa(1:19)//'.INFO'
    call jeveuo(sdappa_info, 'L', vi=v_sdappa_info)
!
! - Get parameters
!
    nb_elem_patch = v_sdappa_info(6*(patch_indx-1)+1)
    do i_elem_patch = 1, nb_elem_patch
        list_elem(i_elem_patch) = v_sdappa_info(6*(patch_indx-1)+1+i_elem_patch)
    end do
!
end subroutine
