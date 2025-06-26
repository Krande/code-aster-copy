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

subroutine endotot(dth1, dflu1, end3d, &
                   wpl3, vwpl33, vwpl33t, wplx3, vwplx33, vwplx33t, &
                   gft, gfr, iso, sigf6, sigf6d, rt33, ref33, &
                   souplesse66, epspg6, eprg00, a, b, x, ipzero, ngf, &
                   ekdc, epspc6, dt3, dr3, dgt3, dgc3, dc, wl3, xmt, dtiso, rt, dtr, &
                   dim3, ndim, ifour, epeqpc, ept, errgf)

!=====================================================================

! Initiation
    implicit none
#include "asterfort/endo3d.h"
! Entier
    integer(kind=8) i
    integer(kind=8) ngf, ndim, ifour
    integer(kind=8) ipzero(ngf)
! Carac
    aster_logical ::  end3d, iso, dtiso
! Reel
    real(kind=8) :: dth1, dflu1
    real(kind=8) :: dtr, dt3(3)
    real(kind=8) :: wpl3(3), vwpl33(3, 3), vwpl33t(3, 3), wplx3(3), vwplx33(3, 3), vwplx33t(3, 3)
    real(kind=8) :: gft, gfr, sigf6(6), sigf6d(6), rt33(3, 3), ref33(3, 3)
    real(kind=8) :: souplesse66(6, 6), epspg6(6), eprg00, a(ngf, ngf+1), b(ngf), x(ngf)
    real(kind=8) :: ekdc, epspc6(6), dr3(3), dgt3(3), dgc3(3), dc, wl3(3), xmt, rt
    real(kind=8) :: dim3, epeqpc, ept, errgf
    real(kind=8) :: umdt

!***********************************************************************
!       prise en compte de l'endommagement mécanique
    if (end3d) then
!            calcul des endommagements et ouvertures de fissures
        call endo3d(wpl3, vwpl33, vwpl33t, wplx3, vwplx33, vwplx33t, &
                    gft, gfr, iso, sigf6, sigf6d, rt33, ref33, &
                    souplesse66, epspg6, eprg00, a, b, x, ipzero, ngf, &
                    ekdc, epspc6, dt3, dr3, dgt3, dgc3, dc, wl3, xmt, dtiso, rt, dtr, &
                    dim3, ndim, ifour, epeqpc, ept, errgf)
    end if

!***********************************************************************
!   prise en compte de l endo thermique et de fluage
    umdt = (1.d0-dth1)*(1.d0-dflu1)
    do i = 1, 6
        sigf6d(i) = (sigf6d(i))*umdt
    end do

!***********************************************************************
end
