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

subroutine select_dof_gene(nume_equa_genez, nb_cmp, cata_cmp, list_cmp, list_equa, &
                           tabl_equa)
!
    implicit none
!
#include "asterfort/nueq_chck.h"
#include "asterfort/jeveuo.h"
#include "asterfort/assert.h"
!
!
    character(len=*), intent(in) :: nume_equa_genez
    integer(kind=8), intent(in) :: nb_cmp
    character(len=8), pointer, optional :: cata_cmp(:)
    character(len=8), pointer, optional :: list_cmp(:)
    integer(kind=8), pointer, optional :: list_equa(:)
    integer(kind=8), pointer, optional :: tabl_equa(:, :)
!
! --------------------------------------------------------------------------------------------------
!
! Select dof from list of nodes and components from NUME_EQUA_GENE
!
! --------------------------------------------------------------------------------------------------
!
! Output:
!    list_equa    : vector on complete numbering [1:nb_equa]
!                   for ieq =  [1:nb_equa]
!                      list_equa[ieq] = 0 if node+component not present
!                      list_equa[ieq] = 1 if node+component is present
!    tabl_equa    : table on complete numbering [1:nb_equa, 1:nb_cmp]
!                   for ieq = [1:nb_equa]
!                      for icmp = [1:nb_cmp]
!                         tabl_equa[ieq,icmp] = 0 if node+component not present
!                         tabl_equa[ieq,icmp] = 1 if node+component is present
!
! In  nume_equa_gene     : name of profile (NUME_EQUA_GENE)
! IO  list_equa     : list of equations
! IO  tabl_equa     : table of equations by components
! In  nb_cmp        : number of components
! In  cata_cmp      : list of components in catalog (name)
! In  list_cmp      : list of components (name)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_equa, i_cmp, nb_equa
    integer(kind=8) :: node_nume
    character(len=8) :: name_cmp
    character(len=19) :: nume_equa_gene
    integer(kind=8), pointer :: v_desc(:) => null()
    integer(kind=8), pointer :: v_deeq(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nume_equa_gene = nume_equa_genez
    call nueq_chck(nume_equa_gene, nb_equa)
    call jeveuo(nume_equa_gene//'.DESC', 'L', vi=v_desc)
    ASSERT(v_desc(1) .eq. 2)
!
    call jeveuo(nume_equa_gene//'.DEEQ', 'L', vi=v_deeq)
    do i_equa = 1, nb_equa
        node_nume = v_deeq(2*i_equa)
        do i_cmp = 1, nb_cmp
            if (present(list_cmp)) then
                name_cmp = list_cmp(i_cmp)
            else
                name_cmp = cata_cmp(i_cmp)
            end if
            if (name_cmp .eq. 'LAGR' .and. node_nume .lt. 0) then
                if (present(tabl_equa)) then
                    tabl_equa(i_equa, i_cmp) = 1
                elseif (present(list_equa)) then
                    list_equa(i_equa) = 1
                else
                    ASSERT(.false.)
                end if
            end if
            if (name_cmp .eq. 'GENE' .and. node_nume .gt. 0) then
                if (present(tabl_equa)) then
                    tabl_equa(i_equa, i_cmp) = 1
                elseif (present(list_equa)) then
                    list_equa(i_equa) = 1
                else
                    ASSERT(.false.)
                end if
            end if
        end do
    end do

end subroutine
