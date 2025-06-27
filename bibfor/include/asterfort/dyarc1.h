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
    subroutine dyarc1(instc, nbpas, insta, nbinst, arch,&
                      epsi, crit)
        real(kind=8) :: instc(*)
        integer(kind=8) :: nbpas
        real(kind=8) :: insta(*)
        integer(kind=8) :: nbinst
        integer(kind=8) :: arch(*)
        real(kind=8) :: epsi
        character(len=8) :: crit
    end subroutine dyarc1
end interface
