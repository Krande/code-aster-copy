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
    subroutine lcmmlc(nmat, nbcomm, cpmono, nfs, nsg,&
                      hsr, nsfv, nsfa, ifa, nbsys,&
                      is, dt, nvi, vind, yd,&
                      dy, itmax, toler, materf, expbp,&
                      taus, dalpha, dgamma, dp, crit,&
                      sgns, rp, iret)
        integer(kind=8) :: nvi
        integer(kind=8) :: nsg
        integer(kind=8) :: nmat
        integer(kind=8) :: nbcomm(nmat, 3)
        character(len=24) :: cpmono(5*nmat+1)
        integer(kind=8) :: nfs
        real(kind=8) :: hsr(nsg, nsg)
        integer(kind=8) :: nsfv
        integer(kind=8) :: nsfa
        integer(kind=8) :: ifa
        integer(kind=8) :: nbsys
        integer(kind=8) :: is
        real(kind=8) :: dt
        real(kind=8) :: vind(nvi)
        real(kind=8) :: yd(*)
        real(kind=8) :: dy(*)
        integer(kind=8) :: itmax
        real(kind=8) :: toler
        real(kind=8) :: materf(nmat*2)
        real(kind=8) :: expbp(nsg)
        real(kind=8) :: taus
        real(kind=8) :: dalpha
        real(kind=8) :: dgamma
        real(kind=8) :: dp
        real(kind=8) :: crit
        real(kind=8) :: sgns
        real(kind=8) :: rp
        integer(kind=8) :: iret
    end subroutine lcmmlc
end interface
