! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
    subroutine nzcomp(jvMaterCode, metaPara, &
                      numeComp, nbPhase, nbVari, &
                      dt10, dt21, inst2, &
                      tno0, tno1, tno2, &
                      metaPrev, metaCurr)
        use Metallurgy_type
        integer, intent(in) :: jvMaterCode
        type(META_MaterialParameters), intent(in) :: metaPara
        integer, intent(in) :: numeComp, nbPhase, nbVari
        real(kind=8), intent(in) :: dt10, dt21, inst2
        real(kind=8), intent(in) :: tno0, tno1, tno2
        real(kind=8), intent(in) :: metaPrev(*)
        real(kind=8), intent(out) :: metaCurr(*)
    end subroutine nzcomp
end interface
