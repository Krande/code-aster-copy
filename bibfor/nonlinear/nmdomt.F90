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
subroutine nmdomt(ds_algopara, ds_algorom_)
!
    use NonLin_Datastructure_type
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterc/getexm.h"
#include "asterfort/assert.h"
#include "asterfort/romAlgoNLRead.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
!
    type(NL_DS_AlgoPara), intent(inout) :: ds_algopara
    type(ROM_DS_AlgoPara), optional, intent(inout) :: ds_algorom_
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm management
!
! Read parameters for algorithm management
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_algopara      : datastructure for algorithm parameters
! IO  ds_algorom       : datastructure for ROM parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: reac_incr, reac_iter, reac_iter_elas
    real(kind=8) :: pas_mini_elas
    integer(kind=8) :: iret
    character(len=16) :: keywf, algo_meth, matrix_pred, matrix_corr, answer
    character(len=8) :: result_prev_disp
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE12_6')
    end if
!
! - Get method
!
    if (getexm(' ', 'METHODE') .eq. 1) then
        call getvtx(' ', 'METHODE', scal=algo_meth)
    else
        algo_meth = 'NEWTON'
    end if
    ds_algopara%method = algo_meth
!
! - Get parameters of method
!
    keywf = 'NEWTON'
    if ((algo_meth .eq. 'NEWTON') .or. (algo_meth .eq. 'NEWTON_KRYLOV')) then
        call getvtx(keywf, 'MATRICE', iocc=1, scal=matrix_corr)
        call getvtx(keywf, 'PREDICTION', iocc=1, scal=matrix_pred, nbret=iret)
        if (iret .eq. 0) then
            matrix_pred = matrix_corr
        end if
        ds_algopara%matrix_pred = matrix_pred
        ds_algopara%matrix_corr = matrix_corr
        if (matrix_pred .eq. 'DEPL_CALCULE') then
            call getvid(keywf, 'EVOL_NOLI', iocc=1, scal=result_prev_disp, nbret=iret)
            if (iret .le. 0) then
                call utmess('F', 'MECANONLINE5_45')
            end if
            ds_algopara%result_prev_disp = result_prev_disp
        end if
        call getvis(keywf, 'REAC_INCR', iocc=1, scal=reac_incr)
        ASSERT(reac_incr .ge. 0)
        ds_algopara%reac_incr = reac_incr
        call getvis(keywf, 'REAC_ITER', iocc=1, scal=reac_iter)
        ASSERT(reac_iter .ge. 0)
        ds_algopara%reac_iter = reac_iter
        call getvr8(keywf, 'PAS_MINI_ELAS', iocc=1, scal=pas_mini_elas, nbret=iret)
        if (iret .ne. 0) then
            ds_algopara%pas_mini_elas = pas_mini_elas
        end if
        call getvis(keywf, 'REAC_ITER_ELAS', iocc=1, scal=reac_iter_elas)
        ASSERT(reac_iter_elas .ge. 0)
        ds_algopara%reac_iter_elas = reac_iter_elas
        call getvtx(keywf, 'MATR_RIGI_SYME', iocc=1, scal=answer)
        ds_algopara%l_matr_rigi_syme = answer .eq. 'OUI'
    else if (algo_meth .eq. 'IMPLEX') then
        ds_algopara%matrix_pred = 'TANGENTE'
        ds_algopara%reac_incr = 1
    else if (algo_meth .eq. 'MODELE_REDUIT') then
        keywf = 'MODELE_REDUIT'
        call getvtx(keywf, 'MATRICE', iocc=1, scal=matrix_corr)
        call getvtx(keywf, 'PREDICTION', iocc=1, scal=matrix_pred, nbret=iret)
        if (iret .eq. 0) then
            matrix_pred = matrix_corr
        end if
        ds_algopara%matrix_pred = matrix_pred
        ds_algopara%matrix_corr = matrix_corr
        if (matrix_pred .eq. 'DEPL_CALCULE') then
            call getvid(keywf, 'EVOL_NOLI', iocc=1, scal=result_prev_disp, nbret=iret)
            if (iret .le. 0) then
                call utmess('F', 'MECANONLINE5_45')
            end if
            ds_algopara%result_prev_disp = result_prev_disp
        end if
        call getvis(keywf, 'REAC_INCR', iocc=1, scal=reac_incr)
        ASSERT(reac_incr .ge. 0)
        ds_algopara%reac_incr = reac_incr
        call getvis(keywf, 'REAC_ITER', iocc=1, scal=reac_iter)
        ASSERT(reac_iter .ge. 0)
        ds_algopara%reac_iter = reac_iter
        call romAlgoNLRead(ds_algorom_)
    else
        ASSERT(.false.)
    end if
!
end subroutine
