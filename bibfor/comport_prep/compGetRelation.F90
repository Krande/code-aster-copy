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
subroutine compGetRelation(factorKeyword, iFactorKeyword, rela_comp)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/deprecated_behavior.h"
#include "asterfort/getvtx.h"
!
    character(len=16), intent(in) :: factorKeyword
    integer(kind=8), intent(in) :: iFactorKeyword
    character(len=16), intent(out) :: rela_comp
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Get type of relation
!
! --------------------------------------------------------------------------------------------------
!
! In  iFactorKeyword   : factor keyword index
! Out rela_comp        : name of behaviour relation
!
! --------------------------------------------------------------------------------------------------
!
    rela_comp = ' '
    call getvtx(factorKeyword, 'RELATION', iocc=iFactorKeyword, scal=rela_comp)
    call deprecated_behavior(rela_comp)
    if ((rela_comp(1:4) .eq. 'META') .and. (rela_comp .ne. 'META_LEMA_ANI')) then
        rela_comp = 'KIT_META'
    end if
!
end subroutine
