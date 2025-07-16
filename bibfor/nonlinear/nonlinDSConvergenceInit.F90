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
subroutine nonlinDSConvergenceInit(ds_conv, sderro, model_)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dismoi.h"
#include "asterfort/GetResi.h"
#include "asterfort/getvr8.h"
#include "asterfort/infdbg.h"
#include "asterfort/jeveuo.h"
#include "asterfort/NonLinear_type.h"
#include "asterfort/SetResi.h"
#include "asterfort/utmess.h"
!
    type(NL_DS_Conv), intent(inout) :: ds_conv
    character(len=24), intent(in) :: sderro
    character(len=24), optional, intent(in) :: model_
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Convergence management
!
! Initializations for convergence management
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_conv          : datastructure for convergence management
! In  sderro           : name of datastructure for events in algorithm
! In  model            : name of model
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=8)  :: exicoq, exipou
    real(kind=8) :: resi_glob_rela
    integer(kind=8) :: iret
    aster_logical :: l_resi_user, l_rela, l_maxi, l_refe, l_comp
    character(len=24) :: eventCONVJv
    integer(kind=8), pointer :: eventCONV(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_7')
    end if

! - No information from user: RESI_GLOB_RELA with 1E-6
    call GetResi(ds_conv, type='RESI_GLOB_RELA', user_para_=resi_glob_rela, &
                 l_resi_test_=l_rela)
    call GetResi(ds_conv, type='RESI_GLOB_MAXI', l_resi_test_=l_maxi)
    call GetResi(ds_conv, type='RESI_REFE_RELA', l_resi_test_=l_refe)
    call GetResi(ds_conv, type='RESI_COMP_RELA', l_resi_test_=l_comp)
    l_resi_user = l_rela .or. l_maxi .or. l_refe .or. l_comp
    if (.not. l_resi_user) then
        call SetResi(ds_conv, type_='RESI_GLOB_RELA', &
                     user_para_=1.d-6, l_resi_test_=ASTER_TRUE)
    end if

! - RESI_REFE_RELA with shell and beam
    if (l_refe .and. present(model_)) then
        call dismoi('EXI_COQUE', model_, 'MODELE', repk=exicoq)
        call dismoi('EXI_POUTRE', model_, 'MODELE', repk=exipou)
        if (exicoq .eq. 'OUI' .and. exipou .eq. 'OUI') then
            call utmess('A', 'MECANONLINE5_32')
        end if
    end if

! - Relaxation of convergence criterion: alarm !
    if (l_rela .and. resi_glob_rela .gt. 1.0001d-4) then
        call utmess('A', 'MECANONLINE5_21')
    end if

! - ARRET=NON: alarm !
    if (.not. ds_conv%l_stop) then
        call utmess('A', 'MECANONLINE5_37')
    end if

! - No NEWTON/PAS_MINI_ELAS parameter => using ITER_GLOB_MAXI instead of ITER_GLOB_ELAS
    call getvr8('NEWTON', 'PAS_MINI_ELAS', iocc=1, nbret=iret)
    if (iret .eq. 0) then
        ds_conv%iter_glob_elas = ds_conv%iter_glob_maxi
    end if

! - For line search
    ds_conv%line_sear_coef = 1.d0
    ds_conv%line_sear_iter = 0

! - For law between residual
    eventCONVJv = sderro(1:19)//'.CONV'
    call jeveuo(eventCONVJv, 'E', vi=eventCONV)
    if (ds_conv%lCritereOr) then
        eventCONV(NB_LOOP+1) = 1
    else
        eventCONV(NB_LOOP+1) = 0
    end if

!
end subroutine
