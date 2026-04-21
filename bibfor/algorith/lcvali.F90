! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine lcvali(materPara, &
                  defoComp, ndim, epsm, deps, &
                  instam, instap, codret)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dscal.h"
!
    type(Material_Para), intent(in) :: materPara
    character(len=16), intent(in) :: defoComp
    integer(kind=8), intent(in) :: ndim
    real(kind=8), intent(in) :: deps(6), epsm(6)
    real(kind=8), intent(in) :: instam, instap
    integer(kind=8), intent(inout) :: codret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbProp = 4
    character(len=16), parameter :: propName(4) = (/'EPSI_MAXI', 'VEPS_MAXI', &
                                                    'TEMP_MINI', 'TEMP_MAXI'/)
    real(kind=8) :: propVale(4)
    integer(kind=8) :: propCode(4)
    integer(kind=8) :: iret1, iret2, iret3, iret
    integer(kind=8) :: ndimsi
    real(kind=8) :: eps(6), epsmax, eps2, vepsm
    real(kind=8) :: veps(6)
    real(kind=8) :: veps2, dt, tmax, tmin, temp
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    iret1 = 0
    iret2 = 0
    iret3 = 0
    ndimsi = 2*ndim
    if (defoComp .ne. 'SIMO_MIEHE') then
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    '+', materPara%jvMaterCode, &
                    ' ', 'VERI_BORNE', &
                    0, ' ', [0.d0], &
                    nbProp, propName, propVale, &
                    propCode, 0)

        if (propCode(1) .eq. 0) then
            epsmax = propVale(1)
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, epsm, b_incx, eps, b_incy)
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 1.d0, deps, b_incx, eps, b_incy)
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            eps2 = sqrt(ddot(b_n, eps, b_incx, eps, b_incy))
            if (eps2 .gt. epsmax) then
                iret1 = 4
            end if
        end if

        if (propCode(2) .eq. 0) then
            vepsm = propVale(2)
            dt = instap-instam
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, deps, b_incx, veps, b_incy)
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            call dscal(b_n, 1.d0/dt, veps, b_incx)
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            veps2 = sqrt(ddot(b_n, veps, b_incx, veps, b_incy))
            if (veps2 .gt. vepsm) then
                iret2 = 4
            end if
        end if

        if (propCode(3) .eq. 0) then
            tmin = propVale(3)
            tmax = propVale(4)
            call rcvarc(' ', 'TEMP', '+', &
                        materPara%schemePara%fami, &
                        materPara%schemePara%kpg, &
                        materPara%schemePara%ksp, &
                        temp, iret)
            if (iret .eq. 0) then
                if ((temp .lt. tmin) .or. (temp .gt. tmax)) then
                    iret3 = 4
                end if
            end if
        end if
        codret = max(iret1, iret2)
        codret = max(codret, iret3)
!
    end if
end subroutine
