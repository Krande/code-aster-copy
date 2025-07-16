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
#include "MeshTypes_type.h"
!
interface
    subroutine ircmva(numcmp, indcmp, ncmpve, ncmprf, nvalec,&
                      nbpg, nbsp, adsv, adsd, adsl,&
                      adsk, cplxFormatZ, tymast, modnum, nuanom,&
                      fieldSupport, val, profas, ideb, ifin,&
                      codret)
        integer(kind=8) :: nbsp
        integer(kind=8) :: nbpg
        integer(kind=8) :: nvalec
        integer(kind=8) :: ncmprf
        integer(kind=8) :: ncmpve
        integer(kind=8) :: numcmp(ncmprf)
        integer(kind=8) :: adsv
        integer(kind=8) :: adsd
        integer(kind=8) :: adsl
        integer(kind=8) :: adsk
        integer(kind=8) :: tymast
        integer(kind=8) :: modnum(MT_NTYMAX)
        integer(kind=8) :: nuanom(MT_NTYMAX, *)
        character(len=8), intent(in) :: fieldSupport
        character(len=*), intent(in) :: cplxFormatZ
        character(len=24), intent(in) :: indcmp
        real(kind=8) :: val(ncmpve, nbsp, nbpg, nvalec)
        integer(kind=8) :: profas(*)
        integer(kind=8) :: ideb
        integer(kind=8) :: ifin
        integer(kind=8) :: codret
    end subroutine ircmva
end interface
