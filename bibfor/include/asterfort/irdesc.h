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
    subroutine irdesc(ifi, nbno, prno, nueq, nec,&
                      dg, ncmpmx, vale, nomcmp, titr,&
                      nomnoe, nomsd, nomsym, ir, numnoe,&
                      lmasu)
        integer(kind=8) :: ifi
        integer(kind=8) :: nbno
        integer(kind=8) :: prno(*)
        integer(kind=8) :: nueq(*)
        integer(kind=8) :: nec
        integer(kind=8) :: dg(*)
        integer(kind=8) :: ncmpmx
        complex(kind=8) :: vale(*)
        character(len=*) :: nomcmp(*)
        character(len=*) :: titr
        character(len=*) :: nomnoe(*)
        character(len=*) :: nomsd
        character(len=*) :: nomsym
        integer(kind=8) :: ir
        integer(kind=8) :: numnoe(*)
        aster_logical :: lmasu
    end subroutine irdesc
end interface
