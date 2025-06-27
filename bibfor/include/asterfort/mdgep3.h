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
    subroutine mdgep3(neq, nbexci, psidel, temps, nomfon,&
                      tab, kacce, kprof, inst, indice, taille)
        integer(kind=8), intent(in) :: nbexci
        integer(kind=8), intent(in) :: neq
        real(kind=8), intent(in) :: psidel(neq, nbexci)
        real(kind=8), intent(in) :: temps
        character(len=8), intent(in) :: nomfon(2*nbexci)
        real(kind=8), intent(out) :: tab(neq)
        character(len=4), intent(in), optional :: kacce
        character(len=4), intent(in), optional :: kprof
        integer(kind=8), intent(in), optional :: inst
        integer(kind=8), intent(inout), optional :: indice
        integer(kind=8), intent(inout), optional :: taille
    end subroutine mdgep3
end interface
