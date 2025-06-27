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
!
interface
!
    subroutine cnscno(cnsz, numeqz, prol0, basez, cnoz, &
                      kstop, iret, nbz, vchamz, lprofconst, prolong)
        !
        use proj_champ_module
        !
        character(len=*) :: cnsz
        character(len=*) :: numeqz
        character(len=*) :: prol0
        character(len=*) :: basez
        character(len=*) :: cnoz
        character(len=1) :: kstop
        integer(kind=8)          :: iret
        integer(kind=8),            optional :: nbz
        character(len=24),  optional :: vchamz
        aster_logical,      optional :: lprofconst
        type(prolongation), optional :: prolong
    end subroutine cnscno
!
end interface
