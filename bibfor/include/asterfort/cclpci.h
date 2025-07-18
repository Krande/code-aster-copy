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
    subroutine cclpci(option, &
                      modelZ, materFieldZ, materCodeZ, caraElemZ, &
                      resultIn, resultOut, &
                      ligrel, numeStore, &
                      nbParaIn, lpain, lchin)
        character(len=16), intent(in) :: option
        character(len=*), intent(in) :: modelZ, materFieldZ, materCodeZ, caraElemZ
        character(len=8), intent(in) :: resultIn, resultOut
        character(len=24), intent(in) :: ligrel
        integer(kind=8), intent(in) :: numeStore
        integer(kind=8), intent(out) :: nbParaIn
        character(len=8), intent(out) :: lpain(100)
        character(len=24), intent(out) :: lchin(100)
    end subroutine cclpci
end interface
