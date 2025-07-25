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
subroutine romMultiParaProdModeInit(ds_multipara, nb_mode_maxi)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcrem.h"
#include "asterfort/wkvect.h"
#include "asterfort/codent.h"
!
    type(ROM_DS_MultiPara), intent(inout) :: ds_multipara
    integer(kind=8), intent(in) :: nb_mode_maxi
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction
!
! Initializations for products matrix x mode, reduced matrix and reduced vector
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_multipara     : datastructure for multiparametric problems
! In  nb_mode_maxi     : maximum number of empirical modes
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: i_matr, nb_matr, i_vect, nb_vect, jv_dummy, nbEqua
    character(len=24) :: prod_matr_mode, matr_mode_curr
    character(len=8) :: matr_name, matr_type
    character(len=1) :: prod_type, syst_type
    aster_logical :: l_coefm_cplx
    character(len=7) :: knume
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM2_34')
    end if
!
! - Get parameters
!
    nb_matr = ds_multipara%nb_matr
    nb_vect = ds_multipara%nb_vect
    syst_type = ds_multipara%syst_type
    nbEqua = ds_multipara%field%nbEqua
!
! - Generate names of datastructures
!
    ASSERT(nb_matr .le. 9999999)
    do i_matr = 1, nb_matr
        call codent(i_matr, 'D0', knume)
        ds_multipara%matr_mode_curr(i_matr) = '&&OP0053.MM.'//knume
        ds_multipara%prod_matr_mode(i_matr) = '&&OP0053.PM.'//knume
        ds_multipara%matr_redu(i_matr) = '&&OP0053.MR.'//knume
    end do
    ASSERT(nb_vect .le. 9999999)
    do i_vect = 1, nb_vect
        call codent(i_vect, 'D0', knume)
        ds_multipara%vect_redu(i_vect) = '&&OP0053.VR.'//knume
    end do
!
! - Prepare product [Matrix] x [Mode]
!
    do i_matr = 1, nb_matr
        matr_mode_curr = ds_multipara%matr_mode_curr(i_matr)
        matr_name = ds_multipara%matr_name(i_matr)
        matr_type = ds_multipara%matr_type(i_matr)
        l_coefm_cplx = ds_multipara%matr_coef(i_matr)%l_cplx
        prod_type = 'R'
        if (matr_type .eq. 'C' .or. l_coefm_cplx .or. syst_type .eq. 'C') then
            prod_type = 'C'
        end if
        call vtcrem(matr_mode_curr, matr_name, 'V', prod_type)
        call wkvect(ds_multipara%matr_redu(i_matr), 'V V '//syst_type, nb_mode_maxi*nb_mode_maxi, &
                    jv_dummy)
        prod_matr_mode = ds_multipara%prod_matr_mode(i_matr)
        call wkvect(prod_matr_mode, 'V V '//syst_type, nb_mode_maxi*nbEqua, jv_dummy)
    end do
!
! - Prepare Reduced Vector
!
    do i_vect = 1, nb_vect
        call wkvect(ds_multipara%vect_redu(i_vect), 'V V '//syst_type, nb_mode_maxi, jv_dummy)
    end do
!
end subroutine
