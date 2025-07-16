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
    subroutine utchdl(cham19, nomma, nomail, nonoeu, nupo,&
                      nusp, ivari, nocmp1, iddl, nogranz)
        character(len=*) :: cham19
        character(len=*) :: nomma
        character(len=*) :: nomail
        character(len=*) :: nonoeu
        integer(kind=8) :: nupo
        integer(kind=8) :: nusp
        integer(kind=8) :: ivari
        character(len=*) :: nocmp1
        integer(kind=8) :: iddl
        aster_logical, intent(in), optional :: nogranz
    end subroutine utchdl
end interface
