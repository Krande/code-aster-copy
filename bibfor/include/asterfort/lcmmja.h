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
    subroutine lcmmja(typmod, nmat, materf, timed,&
                      timef, itmax, toler, nbcomm, cpmono,&
                      pgl, nfs, nsg, toutms, hsr,&
                      nr, nvi, vind, df, yf,&
                      yd, dy, drdy, iret)
        integer(kind=8) :: nr
        integer(kind=8) :: nsg
        integer(kind=8) :: nfs
        integer(kind=8) :: nmat
        character(len=8) :: typmod
        real(kind=8) :: materf(nmat*2)
        real(kind=8) :: timed
        real(kind=8) :: timef
        integer(kind=8) :: itmax
        real(kind=8) :: toler
        integer(kind=8) :: nbcomm(nmat, 3)
        character(len=24) :: cpmono(5*nmat+1)
        real(kind=8) :: pgl(3, 3)
        real(kind=8) :: toutms(nfs, nsg, 6)
        real(kind=8) :: hsr(nsg, nsg)
        integer(kind=8) :: nvi
        real(kind=8) :: vind(*)
        real(kind=8) :: df(3, 3)
        real(kind=8) :: yf(*)
        real(kind=8) :: yd(*)
        real(kind=8) :: dy(*)
        real(kind=8) :: drdy(nr, nr)
        integer(kind=8) :: iret
    end subroutine lcmmja
end interface
