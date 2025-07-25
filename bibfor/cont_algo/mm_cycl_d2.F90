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

subroutine mm_cycl_d2(ds_contact, i_cont_poin, &
                      indi_cont_eval, indi_frot_eval, &
                      coef_frot_prev, &
                      pres_frot_curr, pres_frot_prev, &
                      dist_frot_curr, dist_frot_prev, &
                      alpha_frot_matr, alpha_frot_vect)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/iscode.h"
#include "asterfort/iscycl.h"
#include "asterfort/isdeco.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mm_cycl_zonf.h"
#include "asterfort/mm_cycl_erase.h"
#include "asterfort/mm_cycl_shift.h"
#include "asterfort/mm_cycl_laugf.h"
!
! person_in_charge: ayaovi-dzifa.kudawoo at edf.fr
!
    type(NL_DS_Contact), intent(in) :: ds_contact
    integer(kind=8), intent(in) :: i_cont_poin
    integer(kind=8), intent(in) :: indi_cont_eval
    integer(kind=8), intent(in) :: indi_frot_eval
    real(kind=8), intent(in)  :: pres_frot_curr(3), pres_frot_prev(3)
    real(kind=8), intent(in)  :: dist_frot_curr(3), dist_frot_prev(3)
    real(kind=8), intent(in)  :: coef_frot_prev
    real(kind=8), intent(out)  :: alpha_frot_matr, alpha_frot_vect
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve - Cycling
!
! Detection: sliding/sticking
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_contact       : datastructure for contact management
! In  i_cont_poin      : contact point index
! In  indi_cont_eval   : evaluation of new contact status
! In  indi_frot_eval   : evaluation of new friction status
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: sdcont_cyclis
    integer(kind=8), pointer :: p_sdcont_cyclis(:) => null()
    character(len=24) :: sdcont_cycnbr
    integer(kind=8), pointer :: p_sdcont_cycnbr(:) => null()
    character(len=24) :: sdcont_cyceta
    integer(kind=8), pointer :: p_sdcont_cyceta(:) => null()
    integer(kind=8) :: statut(30)
    integer(kind=8) :: cycl_type, cycl_long_acti
    integer(kind=8) :: cycl_ecod(1), cycl_long, cycl_stat
    aster_logical :: detect
    integer(kind=8)       :: zone_frot_prev, zone_frot_curr
    real(kind=8)  :: nrese_curr, nrese_prev
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    cycl_long_acti = 3
    cycl_type = 2
    detect = .false.
!
! - Access to cycling objects
!
    sdcont_cyclis = ds_contact%sdcont_solv(1:14)//'.CYCLIS'
    sdcont_cycnbr = ds_contact%sdcont_solv(1:14)//'.CYCNBR'
    sdcont_cyceta = ds_contact%sdcont_solv(1:14)//'.CYCETA'
    call jeveuo(sdcont_cyclis, 'E', vi=p_sdcont_cyclis)
    call jeveuo(sdcont_cycnbr, 'E', vi=p_sdcont_cycnbr)
    call jeveuo(sdcont_cyceta, 'E', vi=p_sdcont_cyceta)
!
! - Cycle state
!
    cycl_ecod(1) = p_sdcont_cyclis(4*(i_cont_poin-1)+cycl_type)
    cycl_long = p_sdcont_cycnbr(4*(i_cont_poin-1)+cycl_type)
    call isdeco(cycl_ecod(1), statut, 30)
!
! - No contact: cycling break
!
    if (indi_cont_eval .eq. 0) then
        call mm_cycl_erase(ds_contact, cycl_type, i_cont_poin)
        goto 99
    end if
!
! - New iteration in cycle
!
    cycl_long = cycl_long+1
    statut(cycl_long) = indi_frot_eval
    call iscode(statut, cycl_ecod(1), 30)
!
! - Cycling detection
!
    cycl_stat = 0
    if (cycl_long .eq. cycl_long_acti) then
        detect = iscycl(cycl_ecod(1), cycl_long_acti)
        if (detect) then
            cycl_stat = 10
!
! - Norm of augmented ratios
!
            call mm_cycl_laugf(pres_frot_prev, dist_frot_prev, coef_frot_prev, nrese_prev)
            call mm_cycl_laugf(pres_frot_curr, dist_frot_curr, coef_frot_prev, nrese_curr)
            call mm_cycl_zonf(nrese_prev, 0.98D0, 1.02D0, zone_frot_prev)
            call mm_cycl_zonf(nrese_curr, 0.98D0, 1.02D0, zone_frot_curr)

            !
            ! - Sub-cycling 1 : grazing adherence
            !
            if (((zone_frot_prev .eq. -1) .and. (zone_frot_curr .eq. 1)) .or. &
                ((zone_frot_prev .eq. 1) .and. (zone_frot_curr .eq. -1))) then
                cycl_stat = 11
                if (zone_frot_prev .eq. 3) then
                    alpha_frot_matr = 0.9
                    alpha_frot_vect = 0.9
                else
                    alpha_frot_matr = 0.1
                    alpha_frot_vect = 0.1
                end if

                !
                ! - Sub-cycling 2
                !
            elseif (((zone_frot_prev .eq. -2) .and. (zone_frot_curr .eq. 1)) .or. &
                    ((zone_frot_prev .eq. 1) .and. (zone_frot_curr .eq. -2))) then
                cycl_stat = 12
                if (zone_frot_prev .eq. 4) then
                    alpha_frot_matr = 1.d-1
                    alpha_frot_vect = 1.d-1
                else
                    alpha_frot_matr = 9.d-1
                    alpha_frot_vect = 9.d-1
                end if

                !
                ! - Sub-cycling 3
                !
            elseif (((zone_frot_prev .eq. -1) .and. (zone_frot_curr .eq. 2)) .or. &
                    ((zone_frot_prev .eq. 2) .and. (zone_frot_curr .eq. -1))) then
                cycl_stat = 13
                if (zone_frot_prev .eq. 3) then
                    alpha_frot_matr = 0.9
                    alpha_frot_vect = 0.9
                else
                    alpha_frot_matr = 0.1
                    alpha_frot_vect = 0.1
                end if

                !
                ! - Sub-cycling 4
                !
            elseif (((zone_frot_prev .eq. 2) .and. (zone_frot_curr .eq. -2)) .or. &
                    ((zone_frot_prev .eq. -2) .and. (zone_frot_curr .eq. 2))) then
                cycl_stat = 14
                alpha_frot_matr = 5.d-1
                alpha_frot_vect = 5.d-1

            else
!                call utmess('A',CONTACT5_2)
            end if
        end if
    end if
!
! - End of cycling detection zone: shifting
!
    if (cycl_long .eq. cycl_long_acti) then
        call mm_cycl_shift(cycl_long_acti, cycl_ecod(1), cycl_long)
    end if
!
! - Cycling save
!
    p_sdcont_cyceta(4*(i_cont_poin-1)+cycl_type) = cycl_stat
    p_sdcont_cyclis(4*(i_cont_poin-1)+cycl_type) = cycl_ecod(1)
    p_sdcont_cycnbr(4*(i_cont_poin-1)+cycl_type) = cycl_long
!
99  continue
!
    call jedema()
end subroutine
