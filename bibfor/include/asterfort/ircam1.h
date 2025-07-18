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
    subroutine ircam1(nofimd, nochmd, existc, ncmprf, numpt,&
                      instan, numord, adsd, adsv, adsl,&
                      adsk, partie, indcmp, ncmpve, ntlcmp,&
                      ntncmp, ntucmp, ntproa, nbimpr, caimpi,&
                      caimpk, typech, nomamd, nomtyp, modnum,&
                      nuanom, lfichUniq, nosdfu, codret)
        integer(kind=8) :: nbimpr
        character(len=*) :: nofimd
        character(len=64) :: nochmd
        integer(kind=8) :: existc
        integer(kind=8) :: ncmprf
        integer(kind=8) :: numpt
        real(kind=8) :: instan
        integer(kind=8) :: numord
        integer(kind=8) :: adsd
        integer(kind=8) :: adsv
        integer(kind=8) :: adsl
        integer(kind=8) :: adsk
        character(len=*) :: partie
        integer(kind=8) :: ncmpve
        character(len=24) :: ntlcmp
        character(len=24) :: ntncmp
        character(len=24) :: ntucmp
        character(len=24) :: ntproa
        character(len=24) :: indcmp
        integer(kind=8) :: caimpi(10, nbimpr)
        character(len=*) :: caimpk(3, nbimpr)
        character(len=8) :: typech
        character(len=*) :: nomamd
        character(len=8) :: nomtyp(*)
        integer(kind=8) :: modnum(MT_NTYMAX), nuanom(MT_NTYMAX, *)
        aster_logical :: lfichUniq
        character(len=8) :: nosdfu
        integer(kind=8) :: codret
    end subroutine ircam1
end interface
