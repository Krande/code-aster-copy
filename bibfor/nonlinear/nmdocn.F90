! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine nmdocn(ds_conv)
!
use NonLin_Datastructure_type
!
implicit none
!
#include "asterf_types.h"
#include "asterc/r8nnem.h"
#include "asterc/r8vide.h"
#include "asterc/getexm.h"
#include "asterfort/assert.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infdbg.h"
#include "asterfort/utmess.h"
#include "asterfort/SetResi.h"
#include "asterfort/SetResiRefe.h"
!
type(NL_DS_Conv), intent(inout) :: ds_conv
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Convergence management
!
! Read parameters for convergence management
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_conv          : datastructure for convergence management
!
! --------------------------------------------------------------------------------------------------
!
    integer :: ifm, niv
    character(len=16) :: keywf
    integer :: iret, iret_rela, iret_maxi, iret_refe, iret_comp, para_inte
    real(kind=8) :: para_real
    character(len=8) :: answer
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE12_8')
    endif
!
! - Initializations
!
    iret_refe = 0
    iret_comp = 0
    keywf     = 'CONVERGENCE'
!
! - Get convergence parameters (maximum iterations)
!
    call getvis(keywf, 'ITER_GLOB_MAXI', iocc=1, scal=para_inte)
    ds_conv%iter_glob_maxi = para_inte
    if (getexm(keywf,'ITER_GLOB_ELAS') .eq. 1) then
        call getvis(keywf, 'ITER_GLOB_ELAS', iocc=1, scal=para_inte)
        ds_conv%iter_glob_elas = para_inte
    endif
!
! - Get convergence parameters (residuals)
!
    call getvr8(keywf, 'RESI_GLOB_RELA', iocc=1, scal=para_real, nbret=iret_rela)
    if (iret_rela .eq. 1) then
        call SetResi(ds_conv   , type_ = 'RESI_GLOB_RELA', &
                    user_para_ = para_real, l_resi_test_ = .true._1)
    endif
    call getvr8(keywf, 'RESI_GLOB_MAXI', iocc=1, scal=para_real, nbret=iret_maxi)
    if (iret_maxi .eq. 1) then
        call SetResi(ds_conv   , type_ = 'RESI_GLOB_MAXI', &
                    user_para_ = para_real, l_resi_test_ = .true._1)
    endif
    if (getexm(keywf,'RESI_COMP_RELA') .eq. 1) then
        call getvr8(keywf, 'RESI_COMP_RELA', iocc=1, scal=para_real, nbret=iret_comp)
        if (iret_comp .eq. 1) then
            call SetResi(ds_conv   , type_ = 'RESI_COMP_RELA', &
                        user_para_ = para_real, l_resi_test_ = .true._1)
        endif
    endif
    if (getexm(keywf,'RESI_REFE_RELA') .eq. 1) then
        call getvr8(keywf, 'RESI_REFE_RELA', iocc=1, scal=para_real, nbret=iret_refe)
        if (iret_refe .eq. 1) then
            call SetResi(ds_conv   , type_ = 'RESI_REFE_RELA', &
                        user_para_ = para_real, l_resi_test_ = .true._1)
        endif
    endif
!
! - Reference residuals
!
    if (iret_refe .eq.1 ) then
        call getvr8(keywf, 'SIGM_REFE', iocc=1, scal=para_real, nbret=iret)
        if (iret .eq. 1) then
            call SetResiRefe(ds_conv   , type_ = 'SIGM_REFE', &
                             user_para_ = para_real, l_refe_test_ = .true._1)
        endif
        call getvr8(keywf, 'EPSI_REFE', iocc=1, scal=para_real, nbret=iret)
        if (iret .eq. 1) then
            call SetResiRefe(ds_conv   , type_ = 'EPSI_REFE', &
                             user_para_ = para_real, l_refe_test_ = .true._1)
        endif
        call getvr8(keywf, 'FLUX_THER_REFE', iocc=1, scal=para_real, nbret=iret)
        if (iret .eq. 1) then
            call SetResiRefe(ds_conv   , type_ = 'FLUX_THER_REFE', &
                             user_para_ = para_real, l_refe_test_ = .true._1)
        endif
        call getvr8(keywf, 'FLUX_HYD1_REFE', iocc=1, scal=para_real, nbret=iret)
        if (iret .eq. 1) then
            call SetResiRefe(ds_conv   , type_ = 'FLUX_HYD1_REFE', &
                             user_para_ = para_real, l_refe_test_ = .true._1)
        endif
        call getvr8(keywf, 'FLUX_HYD2_REFE', iocc=1, scal=para_real, nbret=iret)
        if (iret .eq. 1) then
            call SetResiRefe(ds_conv   , type_ = 'FLUX_HYD2_REFE', &
                             user_para_ = para_real, l_refe_test_ = .true._1)
        endif
        call getvr8(keywf, 'VARI_REFE', iocc=1, scal=para_real, nbret=iret)
        if (iret .eq. 1) then
            call SetResiRefe(ds_conv   , type_ = 'VARI_REFE', &
                             user_para_ = para_real, l_refe_test_ = .true._1)
        endif
        call getvr8(keywf, 'EFFORT_REFE', iocc=1, scal=para_real, nbret=iret)
        if (iret .eq. 1) then
            call SetResiRefe(ds_conv   , type_ = 'EFFORT_REFE', &
                             user_para_ = para_real, l_refe_test_ = .true._1)
        endif
        call getvr8(keywf, 'MOMENT_REFE', iocc=1, scal=para_real, nbret=iret)
        if (iret .eq. 1) then
            call SetResiRefe(ds_conv   , type_ = 'MOMENT_REFE', &
                             user_para_ = para_real, l_refe_test_ = .true._1)
        endif
        call getvr8(keywf, 'DEPL_REFE', iocc=1, scal=para_real, nbret=iret)
        if (iret .eq. 1) then
            call SetResiRefe(ds_conv   , type_ = 'DEPL_REFE', &
                             user_para_ = para_real, l_refe_test_ = .true._1)
        endif
        call getvr8(keywf, 'LAGR_REFE', iocc=1, scal=para_real, nbret=iret)
        if (iret .eq. 1) then
            call SetResiRefe(ds_conv   , type_ = 'LAGR_REFE', &
                             user_para_ = para_real, l_refe_test_ = .true._1)
        endif
        call getvr8(keywf, 'PI_REFE', iocc=1, scal=para_real, nbret=iret)
        if (iret .eq. 1) then
            call SetResiRefe(ds_conv   , type_ = 'PI_REFE', &
                             user_para_ = para_real, l_refe_test_ = .true._1)
        endif
    endif
!
! - Forced convergence
!
    if (getexm(keywf,'ARRET') .eq. 1) then
        call getvtx(keywf, 'ARRET', iocc=1, scal=answer, nbret=iret)
        if (iret .gt. 0) then
            ds_conv%l_stop = answer .eq. 'OUI'
        endif
    endif
!
end subroutine
