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

subroutine romMultiCoefRead(ds_multicoef, keywfact, iocc)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/getexm.h"
#include "asterfort/assert.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvc8.h"
#include "asterfort/getvid.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    type(ROM_DS_MultiCoef), intent(inout) :: ds_multicoef
    character(len=16), intent(in) :: keywfact
    integer(kind=8), intent(in) :: iocc
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction
!
! Read coefficients for multiparametric problems
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_multicoef     : datastructure for multiparametric problems - Coefficients
! In  keywfact         : name of factor keyword
! In  iocc             : index of factor keyword
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_cplx, l_real, l_func, l_cste
    integer(kind=8) :: nocc_coef_r, nocc_coef_c, nocc_func_r, nocc_func_c
!
! --------------------------------------------------------------------------------------------------
!
    nocc_func_r = 0
    nocc_func_c = 0
    if (getexm(keywfact, 'FONC_R') .eq. 1) then
        call getvid(keywfact, 'FONC_R', iocc=iocc, nbret=nocc_func_r)
    end if
    if (getexm(keywfact, 'FONC_C') .eq. 1) then
        call getvid(keywfact, 'FONC_C', iocc=iocc, nbret=nocc_func_c)
    end if
!
! - Get type of coefficient
!
    call getvr8(keywfact, 'COEF_R', iocc=iocc, nbret=nocc_coef_r)
    call getvc8(keywfact, 'COEF_C', iocc=iocc, nbret=nocc_coef_c)
    l_cste = nocc_coef_r .ne. 0 .or. nocc_coef_c .ne. 0
    l_func = nocc_func_r .ne. 0 .or. nocc_func_c .ne. 0
    l_cplx = nocc_coef_c .ne. 0 .or. nocc_func_c .ne. 0
    l_real = nocc_coef_r .ne. 0 .or. nocc_func_r .ne. 0
!
! - Read informations
!
    if (l_func) then
        if (l_real) then
            call getvid(keywfact, 'FONC_R', iocc=iocc, &
                        scal=ds_multicoef%func_name)
        elseif (l_cplx) then
            call getvid(keywfact, 'FONC_C', iocc=iocc, &
                        scal=ds_multicoef%func_name)
        else
            ASSERT(.false.)
        end if
    elseif (l_cste) then
        if (l_real) then
            call getvr8(keywfact, 'COEF_R', iocc=iocc, &
                        scal=ds_multicoef%coef_cste_real)
        elseif (l_cplx) then
            call getvc8(keywfact, 'COEF_C', iocc=iocc, &
                        scal=ds_multicoef%coef_cste_cplx)
        else
            ASSERT(.false.)
        end if
    else
        ASSERT(.false.)
    end if
!
! - Save informations
!
    ds_multicoef%l_cste = l_cste
    ds_multicoef%l_func = l_func
    ds_multicoef%l_cplx = l_cplx
    ds_multicoef%l_real = l_real
!
end subroutine
