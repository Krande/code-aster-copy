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
    subroutine enelpg(fami, iadmat, instan, igau, angl_naut,&
                      compor, f, sigma, nbvari,&
                      vari, enelas)
        character(len=4) :: fami
        integer(kind=8) :: iadmat
        real(kind=8) :: instan
        integer(kind=8) :: igau
        real(kind=8) :: angl_naut(3)
        character(len=16) :: compor(*)
        real(kind=8) :: f(3, 3)
        real(kind=8) :: sigma(6)
        integer(kind=8) :: nbvari
        real(kind=8) :: vari(*)
        real(kind=8) :: enelas
    end subroutine enelpg
end interface
