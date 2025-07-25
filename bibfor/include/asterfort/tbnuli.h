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
    subroutine tbnuli(tabin, npacri, lipacr, vi, vr,&
                      vc, vk, lprec, lcrit, nume)
        character(len=*) :: tabin
        integer(kind=8) :: npacri
        character(len=*) :: lipacr(*)
        integer(kind=8) :: vi(*)
        real(kind=8) :: vr(*)
        complex(kind=8) :: vc(*)
        character(len=*) :: vk(*)
        real(kind=8) :: lprec(*)
        character(len=*) :: lcrit(*)
        integer(kind=8) :: nume
    end subroutine tbnuli
end interface
