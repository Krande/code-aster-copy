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

subroutine nmcoun(mesh, list_func_acti, solver, nume_dof_, matr_asse, &
                  iter_newt, time_curr, hval_incr, hval_algo, hval_veasse, &
                  resi_glob_rela, ds_measure, ds_contact, ctccvg)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisl.h"
#include "asterfort/isfonc.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmcofr.h"
#include "asterfort/nmunil.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(in) :: list_func_acti(*)
    character(len=19), intent(in) :: solver
    character(len=*), intent(in) :: nume_dof_
    character(len=19), intent(in) :: matr_asse
    integer(kind=8), intent(in) :: iter_newt
    real(kind=8), intent(in) :: time_curr
    character(len=19), intent(in) :: hval_incr(*)
    character(len=19), intent(in) :: hval_algo(*)
    character(len=19), intent(in) :: hval_veasse(*)
    real(kind=8), intent(in) :: resi_glob_rela
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(NL_DS_Contact), intent(inout) :: ds_contact
    integer(kind=8), intent(out) :: ctccvg
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Discrete methods - Update displacement (solve contact/unilateral conditions)
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  list_func_acti   : list of active functionnalities
! In  solver           : datastructure for solver parameters
! In  nume_dof         : name of numbering object (NUME_DDL)
! In  matr_asse        : matrix
! In  iter_newt        : index of current Newton iteration
! In  time_curr        : current time
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
! In  hval_veasse      : hat-variable for vectors (node fields)
! In  resi_glob_rela   : current value of RESI_GLOB_RELA
! IO  ds_measure       : datastructure for measure and statistics management
! IO  ds_contact       : datastructure for contact management
! Out ctccvg           : output code for contact algorithm
!                        -1 - No solving
!                         0 - OK
!                        +1 - Maximum contact iteration
!                        +2 - Singular contact matrix
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_unil, l_cont_disc
    aster_logical :: l_all_verif
    character(len=19) :: disp_cumu_inst, disp_iter, cncine
    character(len=19) :: disp_curr
    character(len=14) :: nume_dof
!
! --------------------------------------------------------------------------------------------------
!
    ctccvg = -1
    nume_dof = nume_dof_(1:14)
!
! - Active functionnalities
!
    l_unil = isfonc(list_func_acti, 'LIAISON_UNILATER')
    l_cont_disc = isfonc(list_func_acti, 'CONT_DISCRET')
!
! - Get fields
!
    call nmchex(hval_incr, 'VALINC', 'DEPPLU', disp_curr)
    call nmchex(hval_algo, 'SOLALG', 'DDEPLA', disp_iter)
    call nmchex(hval_algo, 'SOLALG', 'DEPDEL', disp_cumu_inst)
    call nmchex(hval_veasse, 'VEASSE', 'CNCINE', cncine)
!
! - Discrete contact
!
    if (l_cont_disc) then
        l_all_verif = cfdisl(ds_contact%sdcont_defi, 'ALL_VERIF')
        if (.not. l_all_verif) then
            call nmcofr(mesh, disp_curr, disp_cumu_inst, disp_iter, solver, &
                        nume_dof, matr_asse, iter_newt, time_curr, resi_glob_rela, &
                        ds_measure, ds_contact, ctccvg)
        else
            ctccvg = 0
        end if
    end if
!
! - Unilateral condition
!
    if (l_unil) then
        call nmunil(mesh, disp_curr, disp_iter, solver, matr_asse, &
                    cncine, iter_newt, time_curr, ds_contact, nume_dof, ctccvg)
    end if
!
! - Yes for computation
!
    ASSERT(ctccvg .ge. 0)
!
end subroutine
