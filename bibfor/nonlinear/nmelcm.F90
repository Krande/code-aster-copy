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
subroutine nmelcm(model, &
                  ds_material, ds_contact, &
                  ds_constitutive, ds_measure, &
                  hval_incr, hval_algo, &
                  matr_elem)
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
#include "asterfort/memare.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmelco_prep.h"
#include "asterfort/nmrinc.h"
#include "asterfort/nmtime.h"
#include "asterfort/nmvcex.h"
#include "asterfort/reajre.h"
#include "asterfort/utmess.h"
!
    character(len=24), intent(in) :: model
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Contact), intent(in) :: ds_contact
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19), intent(in) :: hval_incr(*), hval_algo(*)
    character(len=19), intent(out) :: matr_elem
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Continue/LAC methods - Compute elementary matrix for contact
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  ds_material      : datastructure for material parameters
! In  ds_contact       : datastructure for contact management
! In  ds_constitutive  : datastructure for constitutive laws management
! IO  ds_measure       : datastructure for measure and statistics management
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
! Out matr_elem        : elementary matrix
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
    aster_logical :: l_cont_lac, l_all_verif
    character(len=19) :: disp_prev, vite_prev, acce_prev, vite_curr, varc_prev, varc_curr
    character(len=19) :: disp_cumu_inst
    character(len=19) :: time_prev, time_curr
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)

! - Begin measure
    call nmtime(ds_measure, 'Init', 'Cont_Elem')
    call nmtime(ds_measure, 'Launch', 'Cont_Elem')

! - Get hat variables
    call nmchex(hval_algo, 'SOLALG', 'DEPDEL', disp_cumu_inst)
    call nmchex(hval_incr, 'VALINC', 'DEPMOI', disp_prev)
    call nmchex(hval_incr, 'VALINC', 'VITMOI', vite_prev)
    call nmchex(hval_incr, 'VALINC', 'ACCMOI', acce_prev)
    call nmchex(hval_incr, 'VALINC', 'VITPLU', vite_curr)
    call nmchex(hval_incr, 'VALINC', 'COMMOI', varc_prev)
    call nmchex(hval_incr, 'VALINC', 'COMPLU', varc_curr)
    call nmvcex('INST', varc_prev, time_prev)
    call nmvcex('INST', varc_curr, time_curr)

! - Get contact parameters
    l_cont_lac = cfdisl(ds_contact%sdcont_defi, 'FORMUL_LAC')
    l_all_verif = cfdisl(ds_contact%sdcont_defi, 'ALL_VERIF')

    if (.not. l_all_verif .and. ((.not. l_cont_lac) .or. ds_contact%nb_cont_pair .ne. 0)) then
! ----- Display
        if (niv .ge. 2) then
            call utmess('I', 'CONTACT5_27')
        end if

! ----- Init fields
        lpain = " "
        lchin = " "
        lpaout = " "
        lchout = " "

! ----- Prepare input fields
        call nmelco_prep('MATR', &
                         model, ds_material, ds_contact, &
                         disp_prev, vite_prev, acce_prev, vite_curr, disp_cumu_inst, &
                         nbin, lpain, lchin, &
                         option, time_prev, time_curr, ds_constitutive)

! ----- <LIGREL> for contact elements
        ligrel = ds_contact%ligrel_elem_cont

! ----- Preparation of elementary matrix
        call detrsd('MATR_ELEM', matr_elem)
        call memare('V', matr_elem, model, 'RIGI_MECA')

! ----- Prepare output fields
        lpaout(1) = 'PMATUNS'
        lchout(1) = matr_elem(1:15)//'.M01'
        lpaout(2) = 'PMATUUR'
        lchout(2) = matr_elem(1:15)//'.M02'

! ----- Computation
        call calcul('S', option, ligrel, nbin, lchin, &
                    lpain, nbout, lchout, lpaout, base, &
                    'OUI')

! ----- Copy output fields
        call reajre(matr_elem, lchout(1), base)
        call reajre(matr_elem, lchout(2), base)
    end if

! - End measure
    call nmtime(ds_measure, 'Stop', 'Cont_Elem')
    call nmrinc(ds_measure, 'Cont_Elem')
!
    call jedema()
!
end subroutine
