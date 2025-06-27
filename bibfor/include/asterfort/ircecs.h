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
    subroutine ircecs(ifi, ligrel, nbgrel, longr, ncmpmx,&
                      vale, nomcmp, titr, nomel, loc,&
                      celd, nbnoma, permut, maxnod, typma,&
                      nomsd, nomsym, ir, nbmat, nummai,&
                      lmasu, ncmpu, nucmp)
        integer(kind=8) :: maxnod
        integer(kind=8) :: ifi
        integer(kind=8) :: ligrel(*)
        integer(kind=8) :: nbgrel
        integer(kind=8) :: longr(*)
        integer(kind=8) :: ncmpmx
        complex(kind=8) :: vale(*)
        character(len=*) :: nomcmp(*)
        character(len=*) :: titr
        character(len=*) :: nomel(*)
        character(len=*) :: loc
        integer(kind=8) :: celd(*)
        integer(kind=8) :: nbnoma(*)
        integer(kind=8) :: permut(maxnod, *)
        integer(kind=8) :: typma(*)
        character(len=*) :: nomsd
        character(len=*) :: nomsym
        integer(kind=8) :: ir
        integer(kind=8) :: nbmat
        integer(kind=8) :: nummai(*)
        aster_logical :: lmasu
        integer(kind=8) :: ncmpu
        integer(kind=8) :: nucmp(*)
    end subroutine ircecs
end interface
