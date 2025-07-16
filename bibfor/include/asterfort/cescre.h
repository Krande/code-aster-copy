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
#include "asterf_types.h"
!
interface
    subroutine cescre(basez, cesz, typcez, maz, nomgdz,&
                      ncmpg, licmp, npg, nspt, ncmp, undf0_)
        character(len=*) :: basez
        character(len=*) :: cesz
        character(len=*) :: typcez
        character(len=*) :: maz
        character(len=*) :: nomgdz
        integer(kind=8) :: ncmpg
        character(len=*) :: licmp(*)
        integer(kind=8) :: npg(*)
        integer(kind=8) :: nspt(*)
        integer(kind=8) :: ncmp(*)
        aster_logical, optional, intent(in) :: undf0_
    end subroutine cescre
end interface
