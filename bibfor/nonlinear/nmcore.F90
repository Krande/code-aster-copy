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

subroutine nmcore(sdcrit, sderro, list_func_acti, nume_inst, iter_newt, &
                  line_sear_iter, eta, resi_norm, load_norm, ds_conv)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmcoru.h"
#include "asterfort/nmcrel.h"
#include "asterfort/nmevcv.h"
#include "asterfort/nmlecv.h"
#include "asterfort/nmerge.h"
#include "asterfort/GetResi.h"
#include "asterfort/SetResi.h"
#include "asterfort/nmcore_swap.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=19), intent(in) :: sdcrit
    character(len=24), intent(in) :: sderro
    integer(kind=8), intent(in) :: list_func_acti(*)
    integer(kind=8), intent(in) :: nume_inst
    integer(kind=8), intent(in) :: iter_newt
    integer(kind=8), intent(in) :: line_sear_iter
    real(kind=8), intent(in) :: eta
    real(kind=8), intent(in) :: resi_norm
    real(kind=8), intent(in) :: load_norm
    type(NL_DS_Conv), intent(inout) :: ds_conv
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Convergence management
!
! Evaluate convergence of residuals
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcrit           : name of datastructure to save convergence parameters
! In  sderro           : name of datastructure for error management (events)
! In  list_func_acti   : list of active functionnalities
! In  nume_inst        : index of current time step
! In  iter_newt        : index of current Newton iteration
! In  line_sear_iter   : number of iterations for line search
! In  eta              : coefficient for pilotage (continuation)
! In  resi_norm        : norm of equilibrium residual
! In  load_norm        : norm of exterior loads
! IO  ds_conv          : datastructure for convergence management
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: sdcrit_crtr
    real(kind=8), pointer :: v_sdcrit_crtr(:) => null()
    character(len=16) :: eventType
    real(kind=8) :: load_mini, last_resi_conv, user_para, vale_calc
    integer(kind=8) :: i_resi, nb_resi
    aster_logical :: l_resi_test, l_conv, l_swap_rela_maxi, l_swap_comp_rela
    aster_logical :: cvresi
!
! --------------------------------------------------------------------------------------------------
!
    nb_resi = ds_conv%nb_resi
    cvresi = .true._1
!
! - Get previous convergence informations
!
    sdcrit_crtr = sdcrit(1:19)//'.CRTR'
    call jeveuo(sdcrit_crtr, 'E', vr=v_sdcrit_crtr)
    last_resi_conv = v_sdcrit_crtr(7)
    load_mini = v_sdcrit_crtr(6)
!
! - Event: no convergence
!
    call SetResi(ds_conv, l_conv_=.false._1)
    do i_resi = 1, nb_resi
        eventType = ds_conv%list_resi(i_resi)%eventType
        call nmcrel(sderro, eventType, .false._1)
    end do
!
! - Swap convergence criterias if necessary
!
    call nmcore_swap(sderro, nume_inst, load_norm, load_mini, last_resi_conv, &
                     ds_conv)

! - Check residuals stop criterias
    do i_resi = 1, nb_resi
        vale_calc = ds_conv%list_resi(i_resi)%vale_calc
        user_para = ds_conv%list_resi(i_resi)%user_para
        if (ds_conv%l_resi_test(i_resi)) then
            call nmcoru(vale_calc, user_para, l_conv)
            ! Si on a pas la convergence d'une boucle de point fixe,Pas la peine de vérifier vpene.
            if ((.not. l_conv) .and. (i_resi .eq. 7) .and. (.not. ds_conv%l_stop_pene)) &
                l_conv = .true._1

        else
            l_conv = .true._1
        end if
        ds_conv%list_resi(i_resi)%l_conv = l_conv
    end do
!
! - Save events
!
    do i_resi = 1, nb_resi
        eventType = ds_conv%list_resi(i_resi)%eventType
        l_conv = ds_conv%list_resi(i_resi)%l_conv
        l_resi_test = ds_conv%l_resi_test(i_resi)
        if (l_resi_test) then
            call nmcrel(sderro, eventType,.not. l_conv)
        end if
    end do

! - Event: evaluate convergence of residuals
    call nmevcv(sderro, list_func_acti, 'RESI')
    call nmlecv(sderro, 'RESI', cvresi)
!
! - If swapped: retrieve old convergence system for next step
!
    call nmerge(sderro, 'RESI_MAXR', l_swap_rela_maxi)
    if (l_swap_rela_maxi) then
        call SetResi(ds_conv, type_='RESI_GLOB_RELA', l_resi_test_=.true._1)
        call SetResi(ds_conv, type_='RESI_GLOB_MAXI', l_resi_test_=.false._1)
    end if
    call nmerge(sderro, 'RESI_MAXN', l_swap_comp_rela)
    if (l_swap_comp_rela) then
        call SetResi(ds_conv, type_='RESI_GLOB_RELA', l_resi_test_=.false._1)
        call SetResi(ds_conv, type_='RESI_COMP_RELA', l_resi_test_=.true._1)
    end if
!
! - New minimum exterior load
!
    if ((nume_inst .eq. 1) .and. (iter_newt .eq. 0)) then
        load_mini = load_norm
    else
        if (cvresi .and. (.not. l_swap_rela_maxi)) then
            load_mini = min(load_norm, load_mini)
        end if
    end if
!
! - Save informations
!
    call GetResi(ds_conv, type='RESI_GLOB_RELA', vale_calc_=v_sdcrit_crtr(3))
    call GetResi(ds_conv, type='RESI_GLOB_MAXI', vale_calc_=v_sdcrit_crtr(4))
    call GetResi(ds_conv, type='RESI_REFE_RELA', vale_calc_=v_sdcrit_crtr(8))
    call GetResi(ds_conv, type='RESI_COMP_RELA', vale_calc_=v_sdcrit_crtr(9))
    v_sdcrit_crtr(1) = iter_newt+1
    v_sdcrit_crtr(2) = line_sear_iter
    v_sdcrit_crtr(5) = eta
    v_sdcrit_crtr(6) = load_mini
    if (cvresi) then
        v_sdcrit_crtr(7) = resi_norm
    end if
!
end subroutine
