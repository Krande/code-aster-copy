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
#include "asterf_types.h"
!
interface
    subroutine irchme(ifichi, chanom, partie, nochmd, noresu,&
                      nomsym, typech, numord, nbrcmp, nomcmp,&
                      nbnoec, linoec, nbmaec, limaec, lvarie,&
                      sdcarm, carael, paraListNb, paraListName,&
                      nbCmpDyna, lfichUniq, codret)
        integer(kind=8) :: ifichi
        character(len=19) :: chanom
        character(len=*) :: partie
        character(len=64) :: nochmd
        character(len=8) :: noresu
        character(len=16) :: nomsym
        character(len=8) :: typech
        integer(kind=8) :: numord
        integer(kind=8) :: nbrcmp
        character(len=*) :: nomcmp(*)
        integer(kind=8) :: nbnoec
        integer(kind=8) :: linoec(*)
        integer(kind=8) :: nbmaec
        integer(kind=8) :: limaec(*)
        aster_logical :: lvarie
        character(len=8) :: sdcarm, carael
        integer(kind=8), intent(in) :: paraListNb
        character(len=16), pointer :: paraListName(:)
        integer(kind=8), intent(inout) :: nbCmpDyna
        aster_logical :: lfichUniq
        integer(kind=8) :: codret
    end subroutine irchme
end interface
