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
subroutine romAlgoNLTherResidual(ds_algorom, vec2nd, cnvabt, cnresi, cn2mbr, &
                                 resi_rela, resi_maxi)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
#include "blas/ddot.h"
!
    type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
    character(len=24), intent(in) :: vec2nd, cnvabt, cnresi, cn2mbr
    real(kind=8), intent(out) :: resi_rela, resi_maxi
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction - Solving non-linear problem THERMICS
!
! Evaluate residuals in applying HYPER-REDUCTION
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_algorom       : datastructure for ROM parameters
! In  vec2nd           : applied loads
! In  cnvabt           : BT.T LAMBDA for Dirichlet loads
! In  cnresi           : non-linear residual
! In  cn2mbr           : equilibrium residual (to evaluate convergence)
! Out resi_rela        : value for RESI_GLOB_RELA
! Out resi_maxi        : value for RESI_GLOB_MAXI
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_hrom
    character(len=8) :: resultName
    character(len=19) :: mode
    character(len=24) :: fieldName
    integer(kind=8) :: iEqua, nbEqua, nbMode, iMode, iret
    real(kind=8) :: vnorm, resi
    real(kind=8), pointer :: v_mode(:) => null()
    real(kind=8), pointer :: v_cn2mbr(:) => null()
    real(kind=8), pointer :: v_cn2mbrr(:) => null()
    real(kind=8), pointer :: v_vec2nd(:) => null()
    real(kind=8), pointer :: v_vec2ndr(:) => null()
    real(kind=8), pointer :: v_cnvabt(:) => null()
    real(kind=8), pointer :: v_cnvabtr(:) => null()
    real(kind=8), pointer :: v_cnresi(:) => null()
    real(kind=8), pointer :: v_cnresir(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    resi_rela = 0.d0
    resi_maxi = 0.d0
    vnorm = 0.d0
    resi = 0.d0
!
! - Get parameters
!
    l_hrom = ds_algorom%l_hrom
    resultName = ds_algorom%ds_empi%resultName
    nbEqua = ds_algorom%ds_empi%mode%nbEqua
    nbMode = ds_algorom%ds_empi%nbMode
    fieldName = ds_algorom%ds_empi%mode%fieldName
    ASSERT(ds_algorom%ds_empi%mode%fieldSupp .eq. 'NOEU')
!
! - Access to vectors
!
    call jeveuo(cn2mbr(1:19)//'.VALE', 'E', vr=v_cn2mbr)
    call jeveuo(vec2nd(1:19)//'.VALE', 'L', vr=v_vec2nd)
    call jeveuo(cnvabt(1:19)//'.VALE', 'L', vr=v_cnvabt)
    call jeveuo(cnresi(1:19)//'.VALE', 'L', vr=v_cnresi)
!
! - Create residual
!
    do iEqua = 1, nbEqua
        v_cn2mbr(iEqua) = v_vec2nd(iEqua)-v_cnresi(iEqua)-v_cnvabt(iEqua)
    end do
!
! - Truncation of residual
!
    if (l_hrom) then
        do iEqua = 1, nbEqua
            if (ds_algorom%v_equa_int(iEqua) .eq. 1) then
                v_vec2nd(iEqua) = 0.d0
                v_cnvabt(iEqua) = 0.d0
                v_cnresi(iEqua) = 0.d0
            end if
        end do
    end if
!
! - Product of modes by second member
!
    AS_ALLOCATE(vr=v_cn2mbrr, size=nbMode)
    AS_ALLOCATE(vr=v_vec2ndr, size=nbMode)
    AS_ALLOCATE(vr=v_cnresir, size=nbMode)
    AS_ALLOCATE(vr=v_cnvabtr, size=nbMode)
    do iMode = 1, nbMode
        call rsexch(' ', resultName, fieldName, iMode, mode, &
                    iret)
        call jeveuo(mode(1:19)//'.VALE', 'E', vr=v_mode)
        b_n = to_blas_int(nbEqua)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        v_vec2ndr(iMode) = ddot(b_n, v_mode, b_incx, v_vec2nd, b_incy)
        b_n = to_blas_int(nbEqua)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        v_cnvabtr(iMode) = ddot(b_n, v_mode, b_incx, v_cnvabt, b_incy)
        b_n = to_blas_int(nbEqua)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        v_cnresir(iMode) = ddot(b_n, v_mode, b_incx, v_cnresi, b_incy)
    end do
!
! - Compute maximum
!
    do iMode = 1, nbMode
        v_cn2mbrr(iMode) = v_vec2ndr(iMode)-v_cnresir(iMode)-v_cnvabtr(iMode)
        resi = resi+(v_cn2mbrr(iMode))**2
        vnorm = vnorm+(v_vec2ndr(iMode)-v_cnvabtr(iMode))**2
        resi_maxi = max(resi_maxi, abs(v_cn2mbrr(iMode)))
    end do
!
! - Compute relative
!
    if (vnorm .gt. 0.d0) then
        resi_rela = sqrt(resi/vnorm)
    end if
!
! - Cleaning
!
    AS_DEALLOCATE(vr=v_cn2mbrr)
    AS_DEALLOCATE(vr=v_vec2ndr)
    AS_DEALLOCATE(vr=v_cnresir)
    AS_DEALLOCATE(vr=v_cnvabtr)
!
end subroutine
