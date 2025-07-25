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
subroutine epsett(applic, nbrddl, depl, btild, sgmtd, &
                  epsi, wgt, effint)
    implicit none
!
    integer(kind=8) :: nbrddl, i, k
    character(len=6) :: applic
    real(kind=8) :: btild(4, *), depl(*), epsi(*), sgmtd(*), effinb
    real(kind=8) :: wgt, effint(*)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    if (applic .eq. 'DEFORM') then
!
!     CALCULS DES COMPOSANTES DE DEFORMATIONS TRIDIMENSIONNELLES :
!     EPSXX, EPSYY, EPSXY, EPSXZ (CE SONT LES COMPOSANTES TILDE)
!
        do i = 1, 4
            epsi(i) = 0.d0
            do k = 1, nbrddl
                epsi(i) = epsi(i)+btild(i, k)*depl(k)
            end do
        end do
!
    else if (applic .eq. 'EFFORI') then
!
!     CALCULS DES EFFORTS INTERIEURS
!
        do i = 1, nbrddl
            effinb = 0.d0
            do k = 1, 4
                effinb = effinb+btild(k, i)*sgmtd(k)
            end do
            effint(i) = effint(i)+wgt*effinb
        end do
!
    end if
end subroutine
