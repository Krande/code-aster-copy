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
    subroutine lcjacb(fami, kpg, ksp, rela_comp, mod,&
                      nmat, materf, timed, timef,&
                      yf, deps, itmax, toler, nbcomm,&
                      cpmono, pgl, nfs, nsg, toutms,&
                      hsr, nr, nvi, vind,&
                      vinf, epsd, yd, dy,&
                      crit,&
                      drdy, iret)
        integer(kind=8) :: nvi
        integer(kind=8) :: nr
        integer(kind=8) :: nsg
        integer(kind=8) :: nfs
        integer(kind=8) :: nmat
        character(len=*) :: fami
        integer(kind=8) :: kpg
        integer(kind=8) :: ksp
        character(len=16) :: rela_comp
        character(len=8) :: mod
        real(kind=8) :: materf(nmat, 2)
        real(kind=8) :: timed
        real(kind=8) :: timef
        real(kind=8) :: yf(nr)
        real(kind=8) :: deps(*)
        integer(kind=8) :: itmax
        real(kind=8) :: toler
        integer(kind=8) :: nbcomm(nmat, 3)
        character(len=24) :: cpmono(5*nmat+1)
        real(kind=8) :: pgl(3, 3)
        real(kind=8) :: toutms(nfs, nsg, 6)
        real(kind=8) :: hsr(nsg, nsg)
        real(kind=8) :: vind(*)
        real(kind=8) :: vinf(*)
        real(kind=8) :: epsd(*)
        real(kind=8) :: yd(nr)
        real(kind=8) :: dy(nr)
        real(kind=8) :: crit(*)
        real(kind=8) :: drdy(nr, nr)
        integer(kind=8) :: iret
    end subroutine lcjacb
end interface
