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
interface
#include "asterf_types.h"
    subroutine mlclow(fid, lname, gtype, sdim, ecoo,&
                      swm, nip, ipcoo, wght, giname,&
                      isname, cret)
        med_idt :: fid
        character(len=*) :: lname
        med_int :: gtype
        med_int :: sdim
        real(kind=8) :: ecoo(*)
        med_int :: swm
        med_int :: nip
        real(kind=8) :: ipcoo(*)
        real(kind=8) :: wght(*)
        character(len=*) :: giname
        character(len=*) :: isname
        med_int :: cret
    end subroutine mlclow
end interface
