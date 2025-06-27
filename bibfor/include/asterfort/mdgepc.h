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
    subroutine mdgepc(neq, nbmode, bmodal, xgene, u, kacce, kprof, &
                      inst, instt, indice, taille, kcham)
        integer(kind=8), intent(in) :: nbmode
        integer(kind=8), intent(in) :: neq
        real(kind=8), intent(in) :: bmodal(neq, nbmode)
        complex(kind=8), intent(in) :: xgene(nbmode)
        complex(kind=8), intent(out) :: u(neq)
        character(len=4), intent(in), optional :: kacce
        character(len=4), intent(in), optional :: kprof
        integer(kind=8), intent(in), optional :: inst
        integer(kind=8), intent(in), optional :: instt
        integer(kind=8), intent(inout), optional :: indice
        integer(kind=8), intent(inout), optional :: taille
        character(len=24), intent(in), optional :: kcham
    end subroutine mdgepc
end interface
