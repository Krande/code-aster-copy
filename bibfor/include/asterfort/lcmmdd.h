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
    subroutine lcmmdd(taus, coeft, ifa, nmat, nbcomm,&
                      is, nbsys, nfs, nsg, hsr,&
                      vind, dy, dt, rp, nuecou,&
                      dalpha, dgamma, dp, iret)
        integer(kind=8) :: nsg
        integer(kind=8) :: nmat
        real(kind=8) :: taus
        real(kind=8) :: coeft(nmat)
        integer(kind=8) :: ifa
        integer(kind=8) :: nbcomm(nmat, 3)
        integer(kind=8) :: is
        integer(kind=8) :: nbsys
        integer(kind=8) :: nfs
        real(kind=8) :: hsr(nsg, nsg)
        real(kind=8) :: vind(*)
        real(kind=8) :: dy(12)
        real(kind=8) :: dt
        real(kind=8) :: rp
        integer(kind=8) :: nuecou
        real(kind=8) :: dalpha
        real(kind=8) :: dgamma
        real(kind=8) :: dp
        integer(kind=8) :: iret
    end subroutine lcmmdd
end interface
