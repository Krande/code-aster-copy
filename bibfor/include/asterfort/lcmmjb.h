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
    subroutine lcmmjb(taur, materf, cpmono, ifa, nmat,&
                      nbcomm, dt, nuecou, nsfv, nsfa,&
                      ir, is, nbsys, nfs, nsg,&
                      hsr, vind, dy, iexp, expbp,&
                      itmax, toler, dgsdts, dksdts, dgrdbs,&
                      dkrdbs, iret)
        integer(kind=8) :: nsg
        integer(kind=8) :: nmat
        real(kind=8) :: taur
        real(kind=8) :: materf(nmat*2)
        character(len=24) :: cpmono(5*nmat+1)
        integer(kind=8) :: ifa
        integer(kind=8) :: nbcomm(nmat, 3)
        real(kind=8) :: dt
        integer(kind=8) :: nuecou
        integer(kind=8) :: nsfv
        integer(kind=8) :: nsfa
        integer(kind=8) :: ir
        integer(kind=8) :: is
        integer(kind=8) :: nbsys
        integer(kind=8) :: nfs
        real(kind=8) :: hsr(nsg, nsg)
        real(kind=8) :: vind(*)
        real(kind=8) :: dy(*)
        integer(kind=8) :: iexp
        real(kind=8) :: expbp(nsg)
        integer(kind=8) :: itmax
        real(kind=8) :: toler
        real(kind=8) :: dgsdts
        real(kind=8) :: dksdts
        real(kind=8) :: dgrdbs
        real(kind=8) :: dkrdbs
        integer(kind=8) :: iret
    end subroutine lcmmjb
end interface
