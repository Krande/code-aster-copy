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
subroutine getDensity(jvMater, rho, elasKeywordZ_)
!
    use Behaviour_type
!
    implicit none
!
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/get_elas_id.h"
!
    integer(kind=8), intent(in) :: jvMater
    real(kind=8), intent(out) :: rho
    character(len=*), optional, intent(in) :: elasKeywordZ_
!
! --------------------------------------------------------------------------------------------------
!
! Comportment utility
!
! Get density
!
! --------------------------------------------------------------------------------------------------
!
! In  jvMater          : coded material address
! Out rho              : density
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1), parameter :: poum = '+'
    integer(kind=8), parameter :: kpg = 1, kspg = 1
    real(kind=8) :: valeResu(1)
    integer(kind=8) :: valeIret(1)
    character(len=4), parameter :: inteFami = 'FPG1'
    integer(kind=8) :: elasType
    character(len=16) :: elasKeyword
    real(kind=8), parameter :: r8Dummy = 0.d0
!
! --------------------------------------------------------------------------------------------------
!
    if (present(elasKeywordZ_)) then
        elasKeyword = elasKeywordZ_
    else
        call get_elas_id(jvMater, elasType, elasKeyword)
    end if

    if (elasKeyword .eq. 'ELAS' .or. elasKeyword .eq. 'ELAS_ISTR' &
        .or. elasKeyword .eq. 'ELAS_ORTH') then
        call rcvalb(inteFami, kpg, kspg, &
                    poum, jvMater, ' ', elasKeyword, &
                    0, ' ', [r8Dummy], &
                    1, 'RHO', valeResu, valeIret, 1)
    else
        call utmess('F', 'ELEMENTS_50', sk=elasKeyword)
    end if
    rho = valeResu(1)
!
end subroutine
