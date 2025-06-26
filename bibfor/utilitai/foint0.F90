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
subroutine foint0()
    implicit none
!     REMISE A "ZERO" DU COMMON UTILISE PAR FOINT2
!     ------------------------------------------------------------------
!
    integer(kind=8) :: isvind, isvnxt, svnbpa, svpar, nextsv
    integer(kind=8) :: iaprol, iavale, iapara, luvale
    real(kind=8) :: svresu
    character(len=1) :: svtypf
    character(len=2) :: svprgd
    character(len=24) :: svinte
    character(len=16) :: svnomp
    character(len=19) :: svnomf
    common/ifosav/mxsave, mxpara, svnbpa(4), svpar(10, 4),&
     &                isvnxt, isvind(4), nextsv(4)
    common/jfosav/iaprol(4), iavale(4), iapara(4), luvale(4), lupara(4)
    common/rfosav/svresu(4)
    common/kfosav/svnomp(10, 4), svnomf(4),&
     &                svtypf(4), svprgd(4), svinte(4)
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, lupara, mxpara, mxsave
!-----------------------------------------------------------------------
    mxpara = 10
    mxsave = 4
    isvnxt = mxsave
    do i = 1, mxsave
        svnomf(i) = '????????'
        svresu(i) = 0.d0
        isvind(i) = 1
!
        iaprol(i) = 0
        iavale(i) = 0
        luvale(i) = 0
        iapara(i) = 0
        lupara(i) = 0
    end do
    nextsv(1) = 2
    nextsv(2) = 3
    nextsv(3) = 4
    nextsv(4) = 1
!
end subroutine
