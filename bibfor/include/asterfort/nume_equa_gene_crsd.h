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
    subroutine nume_equa_gene_crsd(nume_equa_genez , base      , nb_equa, nb_sstr, nb_link,&
                             model_genez, gran_namez)
        character(len=*), intent(in) :: nume_equa_genez
        character(len=1), intent(in) :: base
        integer(kind=8), intent(in) :: nb_equa
        integer(kind=8), intent(in) :: nb_sstr
        integer(kind=8), intent(in) :: nb_link
        character(len=*), optional, intent(in) :: model_genez
        character(len=*), optional, intent(in) :: gran_namez
    end subroutine nume_equa_gene_crsd
end interface
