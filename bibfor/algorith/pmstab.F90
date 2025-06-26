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
subroutine pmstab(sigm, sigp, epsm, deps, nbvari, &
                  vim, vip, iforta, instam, instap, &
                  iter, nbpar, nompar, table, vr, &
                  igrad, valimp, imptgt, dsidep, nomvi, &
                  nbvita)
!
! aslint: disable=W1504
    implicit none
!-----------------------------------------------------------------------
!           OPERATEUR    CALC_POINT_MAT STOCKAGE DANS LA TBLE RESULTAT
!-----------------------------------------------------------------------
#include "asterfort/fgequi.h"
#include "asterfort/tbajli.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
    integer(kind=8) :: nbvari, nbpar, i, iter, iforta, igrad, ncmp, imptgt, nbvita
    character(len=4) :: nomeps(6), nomsig(6), nomgrd(9)
    character(len=8) :: k8b, table, vk8(2), nomvi(*)
    character(len=16) :: nompar(*)
    real(kind=8) :: deps(9), sigm(6), sigp(6), epsp(9), epsm(9), epst(9)
    real(kind=8) :: vim(*), vip(*), vr(*), equi(17), valimp(9), sigt(6)
    real(kind=8) :: rac2, instam, instap, dsidep(*)
    complex(kind=8) :: cbid
    blas_int :: b_incx, b_incy, b_n
    data nomeps/'EPXX', 'EPYY', 'EPZZ', 'EPXY', 'EPXZ', 'EPYZ'/
    data nomsig/'SIXX', 'SIYY', 'SIZZ', 'SIXY', 'SIXZ', 'SIYZ'/
    data nomgrd/'F11', 'F12', 'F13', 'F21', 'F22', 'F23', 'F31', 'F32', 'F33'/
!
    cbid = (0.d0, 0.d0)
    rac2 = sqrt(2.d0)
    if (igrad .ne. 0) then
        ncmp = 9
    else
        ncmp = 6
    end if
!
!     STOCKAGE DE LA SOLUTION DANS LA TABLE
    if (ncmp .eq. 6) then
        b_n = to_blas_int(ncmp)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, epsm, b_incx, epsp, b_incy)
        b_n = to_blas_int(ncmp)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, deps, b_incx, epsp, &
                   b_incy)
        b_n = to_blas_int(ncmp)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, epsp, b_incx, epsm, b_incy)
        b_n = to_blas_int(ncmp)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, epsp, b_incx, epst, b_incy)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        call dscal(b_n, 1.d0/rac2, epst(4), b_incx)
    else
        b_n = to_blas_int(ncmp)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, valimp, b_incx, epst, b_incy)
        b_n = to_blas_int(ncmp)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, valimp, b_incx, epsm, b_incy)
    end if
    b_n = to_blas_int(nbvari)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vip, b_incx, vim, b_incy)
    instam = instap
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, sigp, b_incx, sigm, b_incy)
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, sigp, b_incx, sigt, b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    call dscal(b_n, 1.d0/rac2, sigt(4), b_incx)
    call fgequi(sigt, 'SIGM_DIR', 3, equi)
!
    if (iforta .eq. 0) then
!
        b_n = to_blas_int(ncmp)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, epst, b_incx, vr(2), b_incy)
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, sigt, b_incx, vr(ncmp+2), b_incy)
        vr(ncmp+8) = equi(16)
        vr(ncmp+9) = equi(1)
        b_n = to_blas_int(nbvita)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, vip, b_incx, vr(1+ncmp+6+2+1), b_incy)
        vr(1) = instap
        vr(nbpar) = iter
!        ajout KTGT
        if (imptgt .eq. 1) then
            b_n = to_blas_int(36)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, dsidep, b_incx, vr(1+6+6+3+nbvari), b_incy)
        end if
        call tbajli(table, nbpar, nompar, [0], vr, &
                    [cbid], k8b, 0)
!
    else
!
        vr(1) = instap
        vk8(1) = 'EPSI'
        do i = 1, ncmp
            vr(2) = epst(i)
            if (igrad .eq. 0) then
                vk8(2) = nomeps(i)
            else
                vk8(2) = nomgrd(i)
            end if
            call tbajli(table, nbpar, nompar, [0], vr, &
                        [cbid], vk8, 0)
!
        end do
        vk8(1) = 'SIGM'
        do i = 1, 6
            vr(2) = sigt(i)
            vk8(2) = nomsig(i)
            call tbajli(table, nbpar, nompar, [0], vr, &
                        [cbid], vk8, 0)
!
        end do
        vk8(1) = 'SIEQ'
        vr(2) = equi(1)
        vk8(2) = 'VMIS'
        call tbajli(table, nbpar, nompar, [0], vr, &
                    [cbid], vk8, 0)
!
        vr(2) = equi(16)
        vk8(2) = 'TRACE'
        call tbajli(table, nbpar, nompar, [0], vr, &
                    [cbid], vk8, 0)
!
        vk8(1) = 'VARI'
        do i = 1, nbvita
            vr(2) = vip(i)
            vk8(2) = nomvi(i)
!            VK8(2)(1:1)='V'
!            call codent(I,'G',VK8(2)(2:8))
            call tbajli(table, nbpar, nompar, [0], vr, &
                        [cbid], vk8, 0)
        end do
!
    end if
!
!
end subroutine
