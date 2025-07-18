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
subroutine epseff(applic, nb1, depl, btild, sgmtd, &
                  epsi, wgt, effint)
    implicit none
!
    character(len=6) :: applic
    integer(kind=8) :: nb1
    real(kind=8) :: btild(5, *), depl(*), epsi(*), sgmtd(*), effinb(42)
    real(kind=8) :: wgt, effint(*)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, k, nddle
!-----------------------------------------------------------------------
    nddle = 5*nb1+2
    if (applic .eq. 'DEFORM') then
!
!     CALCULS DES COMPOSANTES DE DEFORMATIONS TRIDIMENSIONNELLES :
!     EPSXX, EPSYY, EPSXY, EPSXZ, EPSYZ (CE SONT LES COMPOSANTES TILDE)
!
        do i = 1, 5
            epsi(i) = 0.d0
            do k = 1, nddle
                epsi(i) = epsi(i)+btild(i, k)*depl(k)
            end do
        end do
!
    else if (applic .eq. 'EFFORI') then
!
!     CALCULS DES EFFORTS INTERIEURS
!
        do i = 1, nddle
            effinb(i) = 0.d0
            do k = 1, 5
                effinb(i) = effinb(i)+btild(k, i)*sgmtd(k)
            end do
            effint(i) = effint(i)+wgt*effinb(i)
        end do
!
    end if
end subroutine
