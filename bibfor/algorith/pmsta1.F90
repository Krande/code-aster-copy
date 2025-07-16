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
subroutine pmsta1(sigm, sigp, deps, vim, vip, &
                  nbvari, nbvita, iforta, nbpar, nompar, &
                  vr, igrad, typpar, nomvi, sddisc, &
                  liccvg, itemax, conver, actite)
!
    implicit none
!-----------------------------------------------------------------------
!           OPERATEUR    CALC_POINT_MAT STOCKAGE DANS LA TBLE RESULTAT
!-----------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/detrsd.h"
#include "asterfort/fgequi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/pmevdr.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/wkvect.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
!
    integer(kind=8) :: nbpar, i, nbvari, igrad, ncmp, nbvita, iforta, liccvg(5)
    integer(kind=8) :: actite, jvari
    aster_logical :: itemax, conver
    character(len=4) :: nomeps(6), nomsig(6), nomgrd(9)
    character(len=8) :: k8b, typpar(*), nomvi(*), vk8(2)
    character(len=16) :: nompar(*)
    character(len=19) :: tabinc, sddisc
    complex(kind=8) :: cbid
    real(kind=8) :: deps(9), sigm(6), sigp(6), vim(*), vip(*), vr(*), rac2
    real(kind=8) :: dsig(6)
    real(kind=8) :: depst(9), equi(17)
    blas_int :: b_incx, b_incy, b_n
    data nomeps/'EPXX', 'EPYY', 'EPZZ', 'EPXY', 'EPXZ', 'EPYZ'/
    data nomsig/'SIXX', 'SIYY', 'SIZZ', 'SIXY', 'SIXZ', 'SIYZ'/
    data nomgrd/'F11', 'F12', 'F13', 'F21', 'F22', 'F23', 'F31', 'F32', 'F33'/
!-----------------------------------------------------------------------
!
    call jemarq()
!
    cbid = (0.d0, 0.d0)
    rac2 = sqrt(2.d0)
    if (igrad .ne. 0) then
        ncmp = 9
    else
        ncmp = 6
    end if
!
!     CALCUL DES INCREMENTS POUR NMEVDR
!
    b_n = to_blas_int(ncmp)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, deps, b_incx, depst, b_incy)
    if (igrad .eq. 0) then
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        call dscal(b_n, 1.d0/rac2, depst(4), b_incx)
    end if
!
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, sigp, b_incx, dsig, b_incy)
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, -1.d0, sigm, b_incx, dsig, &
               b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    call dscal(b_n, 1.d0/rac2, dsig, b_incx)
    call fgequi(dsig, 'SIGM_DIR', 3, equi)
!
!
    if (iforta .eq. 0) then
!
        tabinc = '&&OP0033.TABINC'
        call detrsd('TABLE', tabinc)
        call tbcrsd(tabinc, 'V')
        call tbajpa(tabinc, nbpar, nompar, typpar)
!
!        VR CONTIENT L'ACCROISSEMENT DE VARIABLES INTERNES
!        ATTENTION, VR EST LIMITE A  9999 VALEURS
        b_n = to_blas_int(nbvita)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, vip, b_incx, vr(1+ncmp+6+3), b_incy)
        b_n = to_blas_int(nbvita)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, -1.d0, vim, b_incx, vr(1+ncmp+6+3), &
                   b_incy)
!
        b_n = to_blas_int(ncmp)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, depst, b_incx, vr(2), b_incy)
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, dsig, b_incx, vr(ncmp+2), b_incy)
        vr(ncmp+8) = equi(16)
        vr(ncmp+9) = equi(1)
!
        call tbajli(tabinc, nbpar, nompar, [0], vr, &
                    [cbid], k8b, 0)
!
    else
!
        tabinc = '&&OPB033.TABINC'
        call detrsd('TABLE', tabinc)
        call tbcrsd(tabinc, 'V')
        call tbajpa(tabinc, nbpar, nompar, typpar)
!
        call wkvect('&&OP0033.VARI', 'V V R8', nbvita, jvari)
!
!        VR CONTIENT L'ACCROISSEMENT DE VARIABLES INTERNES
        b_n = to_blas_int(nbvita)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, vip, b_incx, zr(jvari), b_incy)
        b_n = to_blas_int(nbvita)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, -1.d0, vim, b_incx, zr(jvari), &
                   b_incy)
!
        vr(1) = 0.d0
        vk8(1) = 'EPSI'
        do i = 1, ncmp
            vr(2) = depst(i)
            if (igrad .eq. 0) then
                vk8(2) = nomeps(i)
            else
                vk8(2) = nomgrd(i)
            end if
            call tbajli(tabinc, nbpar, nompar, [0], vr, &
                        [cbid], vk8, 0)
!
        end do
        vk8(1) = 'SIGM'
        do i = 1, 6
            vr(2) = dsig(i)
            vk8(2) = nomsig(i)
            call tbajli(tabinc, nbpar, nompar, [0], vr, &
                        [cbid], vk8, 0)
!
        end do
        vk8(1) = 'SIEQ'
        vr(2) = equi(1)
        vk8(2) = 'VMIS'
        call tbajli(tabinc, nbpar, nompar, [0], vr, &
                    [cbid], vk8, 0)
!
        vr(2) = equi(16)
        vk8(2) = 'TRACE'
        call tbajli(tabinc, nbpar, nompar, [0], vr, &
                    [cbid], vk8, 0)
!
        vk8(1) = 'VARI'
        do i = 1, nbvita
            vr(2) = zr(jvari-1+i)
            vk8(2) = nomvi(i)
            call tbajli(tabinc, nbpar, nompar, [0], vr, &
                        [cbid], vk8, 0)
        end do
!
    end if
!
!     VERIFICATION DES EVENT-DRIVEN
!
    call pmevdr(sddisc, tabinc, liccvg, itemax, conver, &
                actite)
!
    call jedetr('&&OP0033.VARI')
!
    call jedema()
end subroutine
