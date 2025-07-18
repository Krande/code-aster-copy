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
    subroutine avcipr(nbvec, vectn, vectu, vectv, nbordr,&
                      kwork, sommw, vwork, tdisp, tspaq,&
                      ipgn, nomcri, nomfor, fordef, fatsoc,&
                      proaxe, pseuil, method, ncycl, jvmin,&
                      jvmax, jomin, jomax)
        integer(kind=8) :: tdisp
        integer(kind=8) :: nbordr
        integer(kind=8) :: nbvec
        real(kind=8) :: vectn(3*nbvec)
        real(kind=8) :: vectu(3*nbvec)
        real(kind=8) :: vectv(3*nbvec)
        integer(kind=8) :: kwork
        integer(kind=8) :: sommw
        real(kind=8) :: vwork(tdisp)
        integer(kind=8) :: tspaq
        integer(kind=8) :: ipgn
        character(len=16) :: nomcri
        character(len=16) :: nomfor
        aster_logical :: fordef
        real(kind=8) :: fatsoc
        character(len=16) :: proaxe
        real(kind=8) :: pseuil
        character(len=8) :: method
        integer(kind=8) :: ncycl(nbvec)
        integer(kind=8) :: jvmin
        integer(kind=8) :: jvmax
        integer(kind=8) :: jomin
        integer(kind=8) :: jomax
    end subroutine avcipr
end interface
