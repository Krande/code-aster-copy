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
    subroutine epsett(applic, nbrddl, depl, btild, sgmtd,&
                      epsi, wgt, effint)
        character(len=6) :: applic
        integer(kind=8) :: nbrddl
        real(kind=8) :: depl(*)
        real(kind=8) :: btild(4, *)
        real(kind=8) :: sgmtd(*)
        real(kind=8) :: epsi(*)
        real(kind=8) :: wgt
        real(kind=8) :: effint(*)
    end subroutine epsett
end interface
