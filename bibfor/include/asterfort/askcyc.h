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
    subroutine askcyc(craid, ndim, soumat, beta, ni,&
                      nj, na, axok, liax, nbliax,&
                      libid)
        integer(kind=8) :: nbliax
        complex(kind=8) :: craid(*)
        integer(kind=8) :: ndim
        character(len=24) :: soumat
        real(kind=8) :: beta
        integer(kind=8) :: ni
        integer(kind=8) :: nj
        integer(kind=8) :: na
        aster_logical :: axok
        integer(kind=8) :: liax(nbliax)
        integer(kind=8) :: libid(*)
    end subroutine askcyc
end interface
