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

subroutine mmbclc(mesh, model, iter_newt, nume_inst, ds_measure, &
                  sddisc, sddyna, hval_incr, hval_algo, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cfdisl.h"
#include "asterfort/copisd.h"
#include "asterfort/mmbouc.h"
#include "asterfort/mmchml.h"
#include "asterfort/mmligr.h"
#include "asterfort/mmstat.h"
#include "asterfort/mmctcg.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmrinc.h"
#include "asterfort/nmtime.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    character(len=8), intent(in) :: model
    integer(kind=8), intent(in) :: iter_newt
    integer(kind=8), intent(in) :: nume_inst
    character(len=19), intent(in) :: sddisc
    character(len=19), intent(in) :: sddyna
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19), intent(in) :: hval_incr(*)
    character(len=19), intent(in) :: hval_algo(*)
    type(NL_DS_Contact), intent(inout) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Continue methods - Applying generalized Newton method at Newton's iteration
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! In  iter_newt        : index of current Newton iteration
! In  nume_inst        : index of current time step
! In  sddisc           : datastructure for time discretization
! In  sddyna           : dynamic parameters datastructure
! IO  ds_measure       : datastructure for measure and statistics management
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
! IO  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_all_verif, l_newt_cont, l_newt_geom
    character(len=19) :: sdcont_depgeo
    character(len=19) :: disp_curr, disp_cumu_inst
!
! --------------------------------------------------------------------------------------------------
!
!
! - Get hat variables
!
    call nmchex(hval_incr, 'VALINC', 'DEPPLU', disp_curr)
    call nmchex(hval_algo, 'SOLALG', 'DEPDEL', disp_cumu_inst)
    sdcont_depgeo = ds_contact%sdcont_solv(1:14)//'.DEPG'
!
! - Get contact parameters
!
    l_all_verif = cfdisl(ds_contact%sdcont_defi, 'ALL_VERIF')
    l_newt_cont = cfdisl(ds_contact%sdcont_defi, 'CONT_NEWTON')
    l_newt_geom = cfdisl(ds_contact%sdcont_defi, 'GEOM_NEWTON')
!
! - No computation
!
    if (l_all_verif) then
        call mmbouc(ds_contact, 'Cont', 'Set_Convergence')
    end if
!
! - For generalized Newton method
!
    if (.not. l_all_verif) then
!
! ----- New pairing
!
        if (l_newt_geom) then
            call copisd('CHAMP_GD', 'V', disp_curr, sdcont_depgeo)
            call mmctcg(mesh, ds_contact, ds_measure)
        end if
!
! ----- Start timer for preparation of contact
!
        if (l_newt_cont .or. l_newt_geom) then
            call nmtime(ds_measure, 'Launch', 'Cont_Prep')
        end if
!
! ----- Create contact elements
!
        if (l_newt_geom) then
            call mmligr(mesh, model, ds_contact)
        end if
!
! ----- Management of contact loop
!
        if (l_newt_cont .or. l_newt_geom) then
            call mmstat(mesh, iter_newt, nume_inst, &
                        sddisc, disp_curr, disp_cumu_inst, ds_contact)
        end if
!
! ----- Update input field
!
        if (l_newt_cont .or. l_newt_geom) then
            call mmchml(mesh, ds_contact, sddisc, sddyna, nume_inst)
        end if
!
! ----- Stop timer for preparation of contact
!
        if (l_newt_cont .or. l_newt_geom) then
            call nmtime(ds_measure, 'Stop', 'Cont_Prep')
            call nmrinc(ds_measure, 'Cont_Prep ')
        end if
    end if
!
end subroutine
