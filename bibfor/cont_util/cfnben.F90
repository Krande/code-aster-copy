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

subroutine cfnben(sdcont_defi, enti_indx, enti_type, enti_nb_, enti_jdec_)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=24), intent(in) :: sdcont_defi
    character(len=6), intent(in) :: enti_type
    integer(kind=8), intent(in) :: enti_indx
    integer(kind=8), optional, intent(out) :: enti_nb_
    integer(kind=8), optional, intent(out) :: enti_jdec_
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Utility
!
! Access to connectivity tables
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont_defi      : name of contact definition datastructure (from DEFI_CONTACT)
! In  enti_type        : type of entity
!                        'CONINV' - Inverse connectivity => enti_indx is a node
!                        'CONNEX' - Direct connectivity => enti_indx is an element
! Out enti_nb          : number of entities attached
! Out enti_jdec        : shift in contact datastructure
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: enti_nb, enti_jdec
    character(len=24) :: sdcont_pnomaco
    integer(kind=8), pointer :: p_sdcont_pnomaco(:) => null()
    character(len=24) :: sdcont_pmanoco
    integer(kind=8), pointer :: p_sdcont_pmanoco(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    enti_nb = 0
    enti_jdec = 0
!
! - Access to contact datastructure
!
    sdcont_pnomaco = sdcont_defi(1:16)//'.PNOMACO'
    sdcont_pmanoco = sdcont_defi(1:16)//'.PMANOCO'
    call jeveuo(sdcont_pnomaco, 'L', vi=p_sdcont_pnomaco)
    call jeveuo(sdcont_pmanoco, 'L', vi=p_sdcont_pmanoco)
!
! - Get
!
    if (enti_type .eq. 'CONNEX') then
        enti_nb = p_sdcont_pnomaco(enti_indx+1)-p_sdcont_pnomaco(enti_indx)
        enti_jdec = p_sdcont_pnomaco(enti_indx)
    else if (enti_type .eq. 'CONINV') then
        enti_nb = p_sdcont_pmanoco(enti_indx+1)-p_sdcont_pmanoco(enti_indx)
        enti_jdec = p_sdcont_pmanoco(enti_indx)
    else
        ASSERT(.false.)
    end if
!
    if (present(enti_nb_)) then
        enti_nb_ = enti_nb
    end if
    if (present(enti_jdec_)) then
        enti_jdec_ = enti_jdec
    end if
!
!
end subroutine
