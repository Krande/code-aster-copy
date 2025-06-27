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
!
#include "asterf_types.h"
!
interface
    subroutine ircmpn(nofimd, ncmprf, ncmpve, numcmp, exicmp,&
                      nbvato, nbnoec, linoec, adsl, caimpi,&
                      caimpk, profas, innoce, nosdfu)
        integer(kind=8) :: nbvato
        integer(kind=8) :: ncmprf
        character(len=*) :: nofimd
        integer(kind=8) :: ncmpve
        integer(kind=8) :: numcmp(ncmprf)
        aster_logical :: exicmp(nbvato)
        integer(kind=8) :: nbnoec
        integer(kind=8) :: linoec(*)
        integer(kind=8) :: adsl
        integer(kind=8) :: caimpi(10)
        character(len=80) :: caimpk(3)
        integer(kind=8) :: profas(nbvato)
        integer(kind=8) :: innoce(nbvato)
        character(len=8) :: nosdfu
    end subroutine ircmpn
end interface
