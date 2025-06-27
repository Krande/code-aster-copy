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
    subroutine mmimp4(ifm, noma, nummae, iptm, indcoi,&
                      indcon, indfri, indfrn, lfrot, &
                      lgliss, jeu,  lambdc)
        integer(kind=8) :: ifm
        character(len=8) :: noma
        integer(kind=8) :: nummae
        integer(kind=8) :: iptm
        integer(kind=8) :: indcoi
        integer(kind=8) :: indcon
        integer(kind=8) :: indfri
        integer(kind=8) :: indfrn
        aster_logical :: lfrot
        aster_logical :: lgliss
        real(kind=8) :: jeu
        real(kind=8) :: lambdc
    end subroutine mmimp4
end interface
