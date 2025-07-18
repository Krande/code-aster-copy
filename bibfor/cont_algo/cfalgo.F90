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
subroutine cfalgo(mesh, ds_measure, resi_glob_rela, iter_newt, &
                  solver, nume_dof, matr_asse, disp_iter, &
                  disp_cumu_inst, ds_contact, ctccvg)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/algocg.h"
#include "asterfort/algoco.h"
#include "asterfort/algocp.h"
#include "asterfort/algogl.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfeven.h"
#include "asterfort/cfpost.h"
#include "asterfort/cfprep.h"
#include "asterfort/frogdp.h"
#include "asterfort/infdbg.h"
#include "asterfort/utmess.h"
#include "asterfort/dismoi.h"
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_Measure), intent(inout) :: ds_measure
    real(kind=8), intent(in) :: resi_glob_rela
    integer(kind=8), intent(in) :: iter_newt
    character(len=19), intent(in) :: solver
    character(len=14), intent(in) :: nume_dof
    character(len=19), intent(in) :: matr_asse
    character(len=19), intent(in) :: disp_iter
    character(len=19), intent(in) :: disp_cumu_inst
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
! IO  ds_measure       : datastructure for measure and statistics management
! In  resi_glob_rela   : current value of RESI_GLOB_RELA
! In  iter_newt        : index of current Newton iteration
! In  solver           : datastructure for solver parameters
! In  nume_dof         : name of numbering object (NUME_DDL)
! In  matr_asse        : matrix
! In  disp_iter        : displacement iteration
! In  disp_cumu_inst   : displacement increment from beginning of current time
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
    integer(kind=8) :: algo_cont, algo_frot
    character(len=19) :: syme
    aster_logical :: l_gliss, l_first_geom, l_matr_syme
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT> ... DEBUT DE LA RESOLUTION DU CONTACT'
    end if
!
! - Initializations
!
    ctccvg = 0
!
! - Get contact parameters
!
    algo_cont = cfdisi(ds_contact%sdcont_defi, 'ALGO_CONT')
    algo_frot = cfdisi(ds_contact%sdcont_defi, 'ALGO_FROT')
    l_gliss = cfdisl(ds_contact%sdcont_defi, 'CONT_DISC_GLIS')
!
! - First geometric loop
!
    l_first_geom = ds_contact%l_first_geom
!
! - Preparation of contact solving
!
    call cfprep(mesh, matr_asse, disp_iter, disp_cumu_inst, ds_contact)
!
! - Print
!
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT> ...... DEBUT DE REALISATION DU CALCUL'
    end if
!
! - Type of matrix
!
    call dismoi('TYPE_MATRICE', matr_asse, 'MATR_ASSE', repk=syme)
    l_matr_syme = (syme .eq. 'SYMETRI')
!
! - Select algorithm
!
    if (algo_cont .eq. 4) then
        if (algo_frot .eq. 0) then
            call algocp(ds_measure, ds_contact%sdcont_solv, nume_dof, matr_asse)
        else if (algo_frot .eq. 1) then
            call frogdp(ds_measure, ds_contact%sdcont_solv, nume_dof, matr_asse, resi_glob_rela)
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (algo_cont .eq. 1) then
        if (.not. l_matr_syme) then
            call utmess('F', 'CONTACT_1')
        end if
        if (l_gliss) then
            call algogl(ds_measure, ds_contact%sdcont_defi, ds_contact%sdcont_solv, &
                        solver, matr_asse, mesh, &
                        ctccvg)
        else
            call algoco(ds_measure, ds_contact%sdcont_defi, ds_contact%sdcont_solv, &
                        solver, matr_asse, mesh, &
                        ctccvg)
        end if
    else if (algo_cont .eq. 2) then
        if (.not. l_matr_syme) then
            call utmess('F', 'CONTACT_1')
        end if
        if (algo_frot .eq. 0) then
            call algocg(ds_measure, ds_contact%sdcont_defi, ds_contact%sdcont_solv, &
                        solver, matr_asse, ctccvg)
        else
            ASSERT(ASTER_FALSE)
        end if
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Print
!
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT> ...... FIN DE REALISATION DU CALCUL'
    end if
!
! - Post-treatment
!
    call cfpost(mesh, disp_iter, ds_contact, ctccvg)
!
! - Set events
!
    if (iter_newt .eq. 0) then
        call cfeven('INI', ds_contact)
    end if
    call cfeven('FIN', ds_contact)
!
! - Print
!
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT> ... FIN DE LA RESOLUTION DU CONTACT'
    end if
!
    ASSERT(ctccvg .ge. 0)
!
end subroutine
