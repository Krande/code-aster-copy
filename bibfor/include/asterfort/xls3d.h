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
    subroutine xls3d(callst, grille, jltsv, jltsl, jlnsv,&
                     jlnsl, nbno, jcoor, jcoorg, nbmaf,&
                     jdlima, nbsef, jdlise, jconx1, jconx2,&
                     noma)
        aster_logical :: callst
        aster_logical :: grille
        integer(kind=8) :: jltsv
        integer(kind=8) :: jltsl
        integer(kind=8) :: jlnsv
        integer(kind=8) :: jlnsl
        integer(kind=8) :: nbno
        integer(kind=8) :: jcoor
        integer(kind=8) :: jcoorg
        integer(kind=8) :: nbmaf
        integer(kind=8) :: jdlima
        integer(kind=8) :: nbsef
        integer(kind=8) :: jdlise
        integer(kind=8) :: jconx1
        integer(kind=8) :: jconx2
        character(len=8) :: noma
    end subroutine xls3d
end interface
