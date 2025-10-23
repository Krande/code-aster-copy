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
interface
    subroutine enelpg(fami, jvMaterCode, time, kpg, anglNaut, &
                      relaName, defoComp, &
                      f, sigmEner, &
                      nbVari, vari, &
                      enerElas)
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: jvMaterCode
        real(kind=8), intent(in) :: time, anglNaut(3)
        character(len=16), intent(in) :: relaName, defoComp
        integer(kind=8), intent(in) :: kpg
        real(kind=8), intent(in) :: f(3, 3), sigmEner(6)
        integer(kind=8), intent(in) :: nbVari
        real(kind=8), intent(in) :: vari(nbVari)
        real(kind=8), intent(out) :: enerElas
    end subroutine enelpg
end interface
