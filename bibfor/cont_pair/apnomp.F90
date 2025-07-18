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

subroutine apnomp(sdappa, i_poin, poin_name)
!
    implicit none
!
#include "asterfort/jeveuo.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=19), intent(in) :: sdappa
    integer(kind=8), intent(in) :: i_poin
    character(len=16), intent(out) :: poin_name
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing
!
! Ask datastructure - Name of point
!
! --------------------------------------------------------------------------------------------------
!
! In  sdappa           : name of pairing datastructure
! In  i_poin           : index of point (contact or non-contact)
! Out poin_name        : name of point
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: sdappa_noms
    character(len=16), pointer :: v_sdappa_noms(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    sdappa_noms = sdappa(1:19)//'.NOMS'
    call jeveuo(sdappa_noms, 'L', vk16=v_sdappa_noms)
    poin_name = v_sdappa_noms(i_poin)
!
end subroutine
