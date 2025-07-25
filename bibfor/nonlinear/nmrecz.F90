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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine nmrecz(nume_dof, ds_contact, list_func_acti, &
                  cndiri, cnfint, cnfext, cnsstr, &
                  disp_iter, &
                  func)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmequi.h"
#include "asterfort/isfonc.h"
#include "asterfort/isParallelMesh.h"
!
    integer(kind=8), intent(in) :: list_func_acti(*)
    character(len=24), intent(in) :: nume_dof
    type(NL_DS_Contact), intent(in) :: ds_contact
    character(len=19), intent(in) :: cndiri, cnfint, cnfext, cnsstr, disp_iter
    real(kind=8), intent(out) :: func
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm (line search)
!
! Compute function to minimize for line search
!
! --------------------------------------------------------------------------------------------------
!
! In  nume_dof         : name of numbering object (NUME_DDL)
! In  ds_contact       : datastructure for contact management
! In  list_func_acti   : list of active functionnalities
! In  cndiri           : nodal field for support reaction
! In  cnfint           : nodal field for internal force (integration)
! In  cnfext           : nodal field for external force
! In  cnsstr           : nodal field for sub-structuring force
! In  disp_iter        : displacement iteration
! Out function         : function to minimize for line search
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: cnequi, nommai, nume_equa
    real(kind=8), pointer :: v_cnequi(:) => null()
    integer(kind=8) :: nb_equa, rang
    aster_logical :: l_disp, l_pilo, l_macr, l_parallel_mesh
    real(kind=8), pointer :: v_disp_iter(:) => null()
    integer(kind=8), pointer :: pddl(:) => null()
    real(kind=8), dimension(:), allocatable:: disp, resi
    mpi_int :: mrank
    !
! --------------------------------------------------------------------------------------------------
!
    call dismoi('NB_EQUA', nume_dof, 'NUME_DDL', repi=nb_equa)
    call jeveuo(disp_iter(1:19)//'.VALE', 'L', vr=v_disp_iter)
!
! - hpc context or not
!
    call dismoi('NUME_EQUA', disp_iter, 'CHAM_NO', repk=nume_equa)
    call dismoi('NOM_MAILLA', nume_equa, 'NUME_EQUA', repk=nommai)
    l_parallel_mesh = isParallelMesh(nommai)
    if (l_parallel_mesh) then
        call jeveuo(nume_equa//'.PDDL', 'L', vi=pddl)
        call asmpi_info(rank=mrank)
        rang = to_aster_int(mrank)
    end if
!
! - Compute lack of balance forces
!
    l_disp = ASTER_FALSE
    l_pilo = ASTER_FALSE
    l_macr = isfonc(list_func_acti, 'MACR_ELEM_STAT')
    cnequi = '&&CNCHAR.DONN'
    call nmequi(l_disp, l_pilo, l_macr, cnequi, &
                cnfint, cnfext, cndiri, cnsstr, &
                ds_contact)
    call jeveuo(cnequi(1:19)//'.VALE', 'L', vr=v_cnequi)
!
! - Compute function
!
    if (.not. l_parallel_mesh) then
        func = dot_product(v_disp_iter, v_cnequi)
    else
        disp = pack(v_disp_iter, pddl == rang)
        resi = pack(v_cnequi, pddl == rang)
        func = dot_product(disp, resi)
        call asmpi_comm_vect('MPI_SUM', 'R', scr=func)
        deallocate (disp)
        deallocate (resi)
    end if
!
end subroutine
