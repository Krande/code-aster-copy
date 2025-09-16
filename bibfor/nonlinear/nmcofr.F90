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
subroutine nmcofr(mesh, disp_curr, disp_cumu_inst, disp_iter, solver, &
                  nume_dof, matr_asse, iter_newt, time_curr, resi_glob_rela, &
                  ds_measure, ds_contact, ctccvg)
!
    use NonLin_Datastructure_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cfalgo.h"
#include "asterfort/cfgeom.h"
#include "asterfort/infdbg.h"
#include "asterfort/nmtime.h"
!
    character(len=8), intent(in) :: mesh
    character(len=19), intent(in) :: disp_curr
    character(len=19), intent(in) :: disp_cumu_inst
    character(len=19), intent(in) :: disp_iter
    character(len=19), intent(in) :: solver
    character(len=14), intent(in) :: nume_dof
    character(len=19), intent(in) :: matr_asse
    integer(kind=8), intent(in) :: iter_newt
    real(kind=8), intent(in) :: time_curr
    real(kind=8), intent(in) :: resi_glob_rela
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(NL_DS_Contact), intent(inout) :: ds_contact
    integer(kind=8), intent(out) :: ctccvg
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Discrete methods - Solve contact (pairing and algorithm)
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  disp_curr        : current displacements
! In  disp_iter        : displacement iteration
! In  disp_cumu_inst   : displacement increment from beginning of current time
! In  solver           : datastructure for solver parameters
! In  nume_dof         : name of numbering object (NUME_DDL)
! In  matr_asse        : matrix
! In  iter_newt        : index of current Newton iteration
! In  time_curr        : current time
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
    integer(kind=8) :: ifm, niv
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT> DEBUT DU TRAITEMENT DES CONDITIONS DE CONTACT'
    end if

! - Initializations
    ctccvg = -1

! - Pairing
    call cfgeom(iter_newt, mesh, ds_measure, ds_contact, &
                disp_curr, time_curr)

! - Contact solving
    call nmtime(ds_measure, 'Init', 'Cont_Algo')
    call nmtime(ds_measure, 'Launch', 'Cont_Algo')
    call cfalgo(mesh, ds_measure, resi_glob_rela, &
                solver, nume_dof, matr_asse, disp_iter, &
                disp_cumu_inst, ds_contact, ctccvg)
    call nmtime(ds_measure, 'Stop', 'Cont_Algo')

! - Pairing ended
    ds_contact%l_pair = ASTER_FALSE
    ds_contact%l_first_geom = ASTER_FALSE

! - Yes for computation
    ASSERT(ctccvg .ge. 0)
!
end subroutine
