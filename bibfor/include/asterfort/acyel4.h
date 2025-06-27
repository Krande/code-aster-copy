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
    subroutine acyel4(nmcolz, nomobz, nobl, nobc, okpart,&
                      lilig, nblig, licol, nbcol, cmat,&
                      ndim, ideb, jdeb, beta)
        integer(kind=8) :: ndim
        integer(kind=8) :: nbcol
        integer(kind=8) :: nblig
        character(len=*) :: nmcolz
        character(len=*) :: nomobz
        integer(kind=8) :: nobl
        integer(kind=8) :: nobc
        aster_logical :: okpart
        integer(kind=8) :: lilig(nblig)
        integer(kind=8) :: licol(nbcol)
        complex(kind=8) :: cmat(ndim, ndim)
        integer(kind=8) :: ideb
        integer(kind=8) :: jdeb
        real(kind=8) :: beta
    end subroutine acyel4
end interface
