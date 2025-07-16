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
    subroutine select_dof_gene(nume_equa_genez, nb_cmp, cata_cmp, list_cmp, list_equa,&
                               tabl_equa)
        character(len=*), intent(in) :: nume_equa_genez
        integer(kind=8), intent(in) :: nb_cmp
        character(len=8), pointer, optional :: cata_cmp(:)
        character(len=8), pointer, optional :: list_cmp(:)
        integer(kind=8), pointer, optional :: list_equa(:)
        integer(kind=8), pointer, optional :: tabl_equa(:,:)
    end subroutine select_dof_gene
end interface
