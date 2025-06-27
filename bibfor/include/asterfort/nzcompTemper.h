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
    subroutine nzcompTemper(metaPara, numeComp, &
                            nbVari, nbVariTemper, nbVariPrev, &
                            deltaTime12, &
                            temp1, temp2, &
                            prevMetaIsTemper, &
                            metaPrev, metaCurr, metaCurrTemper)
        use Metallurgy_type
        type(META_MaterialParameters), intent(in) :: metaPara
        integer(kind=8), intent(in) :: numeComp, nbVari, nbVariTemper
        real(kind=8), intent(in) ::  deltaTime12
        real(kind=8), intent(in) :: temp1, temp2
        aster_logical, intent(in) :: prevMetaIsTemper
        real(kind=8), intent(in) :: metaPrev(nbVariPrev)
        real(kind=8), intent(in) :: metaCurr(nbVari)
        real(kind=8), intent(out) :: metaCurrTemper(nbVariTemper)
    end subroutine nzcompTemper
end interface
