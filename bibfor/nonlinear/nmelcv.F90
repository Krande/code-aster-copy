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
subroutine nmelcv(model, &
                  ds_material, ds_contact, ds_constitutive, &
                  disp_prev, vite_prev, &
                  acce_prev, vite_curr, &
                  time_prev, time_curr, &
                  disp_cumu_inst, &
                  vect_elem_cont, vect_elem_fric)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/cfdisl.h"
#include "asterfort/detrsd.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/nmelco_prep.h"
#include "asterfort/vemare.h"
#include "asterfort/reajre.h"
#include "asterfort/utmess.h"
!
    character(len=24), intent(in) :: model
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Contact), intent(in) :: ds_contact
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    character(len=19), intent(in) :: disp_prev, vite_prev, acce_prev, vite_curr
    character(len=19), intent(in) :: time_prev, time_curr
    character(len=19), intent(in) :: disp_cumu_inst
    character(len=19), intent(out) :: vect_elem_cont, vect_elem_fric
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Continue/LAC methods - Compute elementary vectors for contact
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  ds_material      : datastructure for material parameters
! In  ds_contact       : datastructure for contact management
! In  ds_constitutive  : datastructure for constitutive laws management
! In  disp_prev        : displacement at beginning of current time
! In  vite_prev        : speed at beginning of current time
! In  vite_curr        : speed at current time
! In  acce_prev        : acceleration at beginning of current time
! In  time_prev        : previous time
! In  time_curr        : current time
! In  disp_cumu_inst   : displacement increment from beginning of current time
! Out vect_elem_cont   : elementary vectors for contact
! Out vect_elem_cont   : elementary vectors for friction
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8), parameter :: nbout = 2, nbin = 36
    character(len=8) :: lpaout(nbout), lpain(nbin)
    character(len=19) :: lchout(nbout), lchin(nbin)
    character(len=1), parameter :: base = 'V'
    character(len=19) :: ligrel
    character(len=16) :: option
    aster_logical :: l_cont_cont, l_cont_lac
    aster_logical :: l_all_verif
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)

! - Get contact parameters
    l_cont_cont = cfdisl(ds_contact%sdcont_defi, 'FORMUL_CONTINUE')
    l_cont_lac = cfdisl(ds_contact%sdcont_defi, 'FORMUL_LAC')
    l_all_verif = cfdisl(ds_contact%sdcont_defi, 'ALL_VERIF')

! - TYPE DE CONTACT
    if (.not. l_all_verif .and. ((.not. l_cont_lac) .or. ds_contact%nb_cont_pair .ne. 0)) then
! ----- Display
        if (niv .ge. 2) then
            call utmess('I', 'CONTACT5_28')
        end if

! ----- Init fields
        lpain = " "
        lchin = " "
        lpaout = " "
        lchout = " "

! ----- Prepare input fields
        call nmelco_prep('VECT', &
                         model, ds_material, ds_contact, &
                         disp_prev, vite_prev, acce_prev, vite_curr, disp_cumu_inst, &
                         nbin, lpain, lchin, &
                         option, time_prev, time_curr, ds_constitutive)

! ----- <LIGREL> for contact elements
        ligrel = ds_contact%ligrel_elem_cont

! ----- Preparation of elementary vectors
        call detrsd('VECT_ELEM', vect_elem_cont)
        call vemare('V', vect_elem_cont, model)
        call detrsd('VECT_ELEM', vect_elem_fric)
        call vemare('V', vect_elem_fric, model)

! ----- Prepare output fields
        lpaout(1) = 'PVECTCR'
        lchout(1) = vect_elem_cont
        lpaout(2) = 'PVECTFR'
        lchout(2) = vect_elem_fric
! ----- Computation

        call calcul('S', option, ligrel, nbin, lchin, &
                    lpain, nbout, lchout, lpaout, base, &
                    'OUI')

! ----- Copy output fields
        call reajre(vect_elem_cont, lchout(1), base)
        call reajre(vect_elem_fric, lchout(2), base)
    end if
!
    call jedema()
!
end subroutine
