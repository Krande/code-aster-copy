! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

subroutine decode_join(nomjoi, dom1, dom2)
    implicit none
#include "asterf_types.h"
#include "asterf_med.h"
#include "asterfort/assert.h"

!
    character(len=*), intent(in) :: nomjoi
    integer, intent(out) :: dom1, dom2
!
!     ------------------------------------------------------------------
!     DECODAGE DU NOM DU JOIN: 'DOM1  DOM2'
!     ------------------------------------------------------------------
!
    integer :: length, ich, start
    character(len=MED_NAME_SIZE) :: nom
!
    nom = nomjoi
!
    length = 15000
    do ich = 1, MED_NAME_SIZE
        if( nom(ich:ich) .eq. ' ' ) then
            length = ich - 1
            exit
        end if
    end do
    ASSERT(length <= MED_NAME_SIZE)
!
    read(nom(1:length), *) dom1
    nom = nomjoi
!
    start = -1
    do ich = length + 1, MED_NAME_SIZE
        if( nom(ich:ich) .ne. ' ' ) then
            start = ich
            exit
        end if
    end do
    ASSERT(start > length)
    do ich = start, MED_NAME_SIZE
        if( nom(ich:ich) == ' ' ) then
            length = ich - 1
            exit
        end if
    end do
    ASSERT(start <= length)

    read(nom(start:length), *) dom2
!
end subroutine
