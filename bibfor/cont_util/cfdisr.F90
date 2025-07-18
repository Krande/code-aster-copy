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
! person_in_charge: mickael.abbas at edf.fr
!
function cfdisr(sdcont_defi_, question_)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mminfr.h"
!
    real(kind=8) :: cfdisr
    character(len=*), intent(in) :: sdcont_defi_
    character(len=*), intent(in) :: question_
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Utility
!
! Get parameter (real)
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont_defi      : name of contact definition datastructure (from DEFI_CONTACT)
! In  question         : question to select parameter
! Out cfdisr           : value for selected parameter
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: sdcont_defi, question
    character(len=24) :: sdcont_paracr
    real(kind=8), pointer :: v_sdcont_paracr(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    sdcont_defi = sdcont_defi_
    question = question_
    cfdisr = 0.d0
!
! - Access to contact datastructure
!
    sdcont_paracr = sdcont_defi(1:8)//'.PARACR'
    call jeveuo(sdcont_paracr, 'L', vr=v_sdcont_paracr)
!
! - Get parameter
!
    if (question .eq. 'RESI_GEOM') then
        cfdisr = v_sdcont_paracr(1)
    else if (question .eq. 'RESI_FROT') then
        cfdisr = v_sdcont_paracr(2)
    else if (question .eq. 'RESI_CONT') then
        cfdisr = v_sdcont_paracr(7)
    else if (question .eq. 'RESI_ABSO') then
        cfdisr = v_sdcont_paracr(4)
    else if (question .eq. 'COEF_RESI') then
        cfdisr = v_sdcont_paracr(5)
    else if (question .eq. 'ALARME_JEU') then
        cfdisr = mminfr(sdcont_defi, 'ALARME_JEU')
    else if (question .eq. 'PROJ_NEWT_RESI') then
        cfdisr = 1d-4
    else if (question .eq. 'PENE_MAXI') then
        cfdisr = v_sdcont_paracr(6)
        if (cfdisr .le. 0.d0) then
            cfdisr = 1.d-2
        end if
    else
        write (6, *) 'QUESTION: ', question
        ASSERT(.false.)
    end if
!
end function
