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
    subroutine lcmmfe(taus, coeft, materf, ifa, nmat,&
                      nbcomm, necoul, is, nbsys, vind,&
                      dy, rp, alphap, gammap, dt,&
                      dalpha, dgamma, dp, crit, sgns,&
                      nfs, nsg, hsr, iret)
        integer(kind=8) :: nsg
        integer(kind=8) :: nmat
        real(kind=8) :: taus
        real(kind=8) :: coeft(nmat)
        real(kind=8) :: materf(nmat)
        integer(kind=8) :: ifa
        integer(kind=8) :: nbcomm(nmat, 3)
        character(len=16) :: necoul
        integer(kind=8) :: is
        integer(kind=8) :: nbsys
        real(kind=8) :: vind(*)
        real(kind=8) :: dy(*)
        real(kind=8) :: rp
        real(kind=8) :: alphap
        real(kind=8) :: gammap
        real(kind=8) :: dt
        real(kind=8) :: dalpha
        real(kind=8) :: dgamma
        real(kind=8) :: dp
        real(kind=8) :: crit
        real(kind=8) :: sgns
        integer(kind=8) :: nfs
        real(kind=8) :: hsr(nsg, nsg)
        integer(kind=8) :: iret
    end subroutine lcmmfe
end interface
