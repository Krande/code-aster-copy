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
subroutine cont_init(mesh, ds_contact, nume_inst, ds_measure, &
                     sddyna, hval_incr, sdnume, list_func_acti)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterfort/cfdisl.h"
#include "asterfort/isfonc.h"
#include "asterfort/mminit.h"
#include "asterfort/mminit_lac.h"
#include "asterfort/cfinit.h"
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(inout) :: ds_contact
    type(NL_DS_Measure), intent(inout) :: ds_measure
    integer(kind=8), intent(in) :: nume_inst
    character(len=19), intent(in) :: hval_incr(*)
    character(len=19), intent(in) :: sddyna
    integer(kind=8), intent(in) :: list_func_acti(*)
    character(len=19), intent(in) :: sdnume
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! All methods - Initializations for current time step
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! IO  ds_contact       : datastructure for contact management
! In  nume_inst        : index of current step time
! In  hval_incr        : hat-variable for incremental values fields
! IO  ds_measure       : datastructure for measure and statistics management
! In  sddyna           : datastructure for dynamic
! In  sdnume           : name of dof positions datastructure
! In  list_func_acti   : list of active functionnalities
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_cont_disc, l_cont_allv, l_cont_cont, l_cont_lac
!
! --------------------------------------------------------------------------------------------------
!
    l_cont_disc = isfonc(list_func_acti, 'CONT_DISCRET')
    l_cont_cont = isfonc(list_func_acti, 'CONT_CONTINU')
    l_cont_lac = isfonc(list_func_acti, 'CONT_LAC')
    l_cont_allv = cfdisl(ds_contact%sdcont_defi, 'ALL_VERIF')
!
    if (.not. l_cont_allv) then
! ----- For discrete contact
        if (l_cont_disc) then
            call cfinit(ds_contact, nume_inst)
        end if
! ----- For continue contact
        if (l_cont_cont) then
            call mminit(mesh, ds_contact, sddyna, hval_incr, ds_measure, &
                        sdnume, nume_inst, list_func_acti)
        end if
! ----- For continue contact (LAC method)
        if (l_cont_lac) then
            call mminit_lac(mesh, ds_contact, hval_incr, ds_measure, &
                            sdnume, nume_inst)
        end if
    end if
!
end subroutine
