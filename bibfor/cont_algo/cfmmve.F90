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
subroutine cfmmve(mesh, ds_contact, hval_incr, time_curr)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/apcalc.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfmmvc.h"
#include "asterfort/cfmmvs.h"
#include "asterfort/cfpoin.h"
#include "asterfort/cfsans.h"
#include "asterfort/cfveri.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/mmpoin.h"
#include "asterfort/mmveri.h"
#include "asterfort/nmchex.h"
#include "asterfort/mreacg.h"
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(inout) :: ds_contact
    character(len=19), intent(in) :: hval_incr(*)
    real(kind=8), intent(in) :: time_curr
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Post-treatment
!
! Post-treatment for no computation methods
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! IO  ds_contact       : datastructure for contact management
! In  hval_incr        : hat-variable for incremental values fields
! In  time_curr        : current time
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: disp_curr
    aster_logical :: l_cont_cont, l_cont_disc, l_cont_allv
    real(kind=8), pointer :: v_ncomp_jeux(:) => null()
    integer(kind=8), pointer :: v_ncomp_loca(:) => null()
    character(len=16), pointer :: v_ncomp_enti(:) => null()
    integer(kind=8), pointer :: v_ncomp_zone(:) => null()
    integer(kind=8) :: nt_ncomp_poin
!
! --------------------------------------------------------------------------------------------------
!
    l_cont_cont = cfdisl(ds_contact%sdcont_defi, 'FORMUL_CONTINUE')
    l_cont_disc = cfdisl(ds_contact%sdcont_defi, 'FORMUL_DISCRETE')
    l_cont_allv = cfdisl(ds_contact%sdcont_defi, 'ALL_VERIF')
!
! - Get fields
!
    call nmchex(hval_incr, 'VALINC', 'DEPPLU', disp_curr)
!
! - Geometry update
!
    if (l_cont_allv) then
        call mreacg(mesh, ds_contact, disp_curr)
    end if
!
! - Create pairing datastructure
!
    if (l_cont_allv) then
        if (l_cont_cont) then
            call mmpoin(mesh, ds_contact)
        else if (l_cont_disc) then
            call cfpoin(mesh, ds_contact)
        else
            ASSERT(.false.)
        end if
        call apcalc('N_To_S', mesh, ds_contact)
    end if
!
! - Prepare datastructures
!
    call cfmmvc(ds_contact, v_ncomp_jeux, v_ncomp_loca, v_ncomp_enti, v_ncomp_zone, &
                nt_ncomp_poin)
!
! - Evaluate
!
    if (l_cont_cont) then
        call mmveri(mesh, ds_contact, time_curr, nt_ncomp_poin, &
                    v_ncomp_jeux, v_ncomp_loca, v_ncomp_enti, v_ncomp_zone)
    else if (l_cont_disc) then
        call cfveri(mesh, ds_contact, time_curr, nt_ncomp_poin, &
                    v_ncomp_jeux, v_ncomp_loca, v_ncomp_enti, v_ncomp_zone)
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - To print interpenetrations
!
    call cfsans(ds_contact, nt_ncomp_poin, v_ncomp_jeux, v_ncomp_enti, v_ncomp_zone)
!
! - Fill CONT_NOEU datastructure
!
    call cfmmvs(ds_contact, nt_ncomp_poin, v_ncomp_jeux, v_ncomp_loca, v_ncomp_zone)
!
! - Clean
!
    AS_DEALLOCATE(vr=v_ncomp_jeux)
    AS_DEALLOCATE(vi=v_ncomp_loca)
    AS_DEALLOCATE(vk16=v_ncomp_enti)
    AS_DEALLOCATE(vi=v_ncomp_zone)
!
end subroutine
