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
subroutine mppsta(h, ldh, v, ldv, ddlsta, &
                  n, vectt, ddlexc, indico, proj)
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
! METHODE DES PUISSANCES AVEC PROJECTION. OPTIMISATION
! SOUS CONTRAINTES INEGALITES
!
! ----------------------------------------------------------------------
!
! IN  H        : MATRICE REDUITE
! IN  LDH      : NOMBRE DE COEFFICIENTS DE H
! IN  V        : MATRICE DE CHANGEMENT DE BASE
! IN  LDV      : NOMBRE DE COEFFICIENTS DE V
! IN  DDLSTA   : POSITION DES DDL_STAB
! IN  N        : DIMENSION ESPACE GLOBAL
! IN/OUT VECTT : MODE CONVERGE
! IN  DDLEXC   : POSITION DDL IMPOSES
! IN  INDICO   : SI =1 ALORS VECTEUR INITIAL =VECTT
!                SINON RANDOM
! IN  PROJ     : =1 PROJECTION
!                =0 PAS DE PROJECTION DANS LA METHODE
!
! ----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
#include "jeveux.h"
#include "asterfort/jedetr.h"
#include "asterfort/r8inir.h"
#include "asterfort/wkvect.h"
#include "blas/dgemv.h"
#include "blas/dlarnv.h"
#include "blas/dnrm2.h"
#include "blas/dscal.h"
    integer(kind=8) :: n, ldh, ldv, indico, proj
    integer(kind=8) :: ddlsta(n), ddlexc(n)
    real(kind=8) :: h(ldh, ldh), v(ldv, ldh)
    real(kind=8) :: vectt(ldv)
    integer(kind=8) :: i, j, nitmax, num, npro
    integer(kind=4) :: iseed(4)
    integer(kind=8) :: x0, bounds, x
    real(kind=8) :: norm
    real(kind=8) :: temp, epsil, one, zero
    real(kind=8) :: crit2, epsf
    blas_int :: b_incx, b_incy, b_lda, b_m, b_n
    blas_int :: b_idist
    parameter(epsil=1.d-15)
    parameter(crit2=1.d-12)
    parameter(nitmax=40000)
    parameter(one=1.0d+0, zero=0.0d+0)
!
    call wkvect('&&MPPSTA.VECT.TEM1', 'V V R', ldh, x0)
    call wkvect('&&MPPSTA.VECT.TEM2', 'V V R', ldh, bounds)
    call wkvect('&&MPPSTA.VECT.TEM3', 'V V R', ldh, x)
!
    call r8inir(ldh, 0.d0, zr(x0), 1)
    call r8inir(ldh, 0.d0, zr(bounds), 1)
    call r8inir(ldh, 0.d0, zr(x), 1)
!
!     INITIALISATION PAR UN VECT RANDOMDE TAILLE LDH
!
    if (indico .eq. 1) then
        b_lda = to_blas_int(ldv)
        b_m = to_blas_int(ldv)
        b_n = to_blas_int(ldh)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('T', b_m, b_n, one, v, &
                   b_lda, vectt, b_incx, zero, zr(x0), &
                   b_incy)
        b_n = to_blas_int(ldh)
        b_incx = to_blas_int(1)
        temp = dnrm2(b_n, zr(x0), b_incx)
        b_n = to_blas_int(ldh)
        b_incx = to_blas_int(1)
        call dscal(b_n, one/temp, zr(x0), b_incx)
        epsf = crit2
        goto 100
    else
        epsf = epsil
    end if
!
    iseed(1) = 1
    iseed(2) = 3
    iseed(3) = 5
    iseed(4) = 7
!
    b_idist = to_blas_int(2)
    b_n = to_blas_int(ldh)
    call dlarnv(b_idist, iseed, b_n, zr(x0))
!
!     PROJECTION DANS LA BASE INITIALE
!
    b_lda = to_blas_int(ldv)
    b_m = to_blas_int(ldv)
    b_n = to_blas_int(ldh)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dgemv('N', b_m, b_n, one, v, &
               b_lda, zr(x0), b_incx, zero, vectt, &
               b_incy)
!
!     PROJECTION DANS SUR LES DDLS_STAB POSITIFS ET CL
!
    do i = 1, n
        if (ddlsta(i) .eq. 0 .and. proj .eq. 1) then
            if (vectt(i) .lt. zero) then
                vectt(i) = zero
            end if
        else
            vectt(i) = vectt(i)*ddlexc(i)
        end if
    end do
!
!     RETOUR DANS LA BASE D'ARNOLDI
!
    b_lda = to_blas_int(ldv)
    b_m = to_blas_int(ldv)
    b_n = to_blas_int(ldh)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dgemv('T', b_m, b_n, one, v, &
               b_lda, vectt, b_incx, zero, zr(x0), &
               b_incy)
!
!     ON NORME X0
!
    b_n = to_blas_int(ldh)
    b_incx = to_blas_int(1)
    temp = dnrm2(b_n, zr(x0), b_incx)
    b_n = to_blas_int(ldh)
    b_incx = to_blas_int(1)
    call dscal(b_n, one/temp, zr(x0), b_incx)
!
!     ON APPLIQUE LA METHODE DES PUISSANCES
!
100 continue
!
    norm = 1.d0
    do j = 1, nitmax
        npro = 0
        num = 0
        b_lda = to_blas_int(ldh)
        b_m = to_blas_int(ldh)
        b_n = to_blas_int(ldh)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('N', b_m, b_n, one, h, &
                   b_lda, zr(x0), b_incx, zero, zr(x), &
                   b_incy)
        b_lda = to_blas_int(ldv)
        b_m = to_blas_int(ldv)
        b_n = to_blas_int(ldh)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('N', b_m, b_n, one, v, &
                   b_lda, zr(x), b_incx, zero, vectt, &
                   b_incy)
        do i = 1, n
            if (ddlsta(i) .eq. 0 .and. proj .eq. 1) then
                num = num+1
                if (vectt(i) .le. zero) then
                    vectt(i) = zero
                    npro = npro+1
                end if
            else
                vectt(i) = vectt(i)*ddlexc(i)
            end if
        end do
        b_lda = to_blas_int(ldv)
        b_m = to_blas_int(ldv)
        b_n = to_blas_int(ldh)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('T', b_m, b_n, one, v, &
                   b_lda, vectt, b_incx, zero, zr(x), &
                   b_incy)
        b_n = to_blas_int(ldh)
        b_incx = to_blas_int(1)
        temp = dnrm2(b_n, zr(x), b_incx)
        b_n = to_blas_int(ldh)
        b_incx = to_blas_int(1)
        call dscal(b_n, one/temp, zr(x), b_incx)
        do i = 1, ldh
            zr(bounds+i-1) = zr(x0+i-1)-zr(x+i-1)
        end do
        b_n = to_blas_int(ldh)
        b_incx = to_blas_int(1)
        norm = dnrm2(b_n, zr(bounds), b_incx)
!
        do i = 1, ldh
            zr(x0+i-1) = zr(x+i-1)
        end do
!
        if (norm .le. epsf) then
            goto 900
        else if (j .eq. nitmax) then
            goto 900
        end if
    end do
!
900 continue
!
!
!     MENAGE
!
    call jedetr('&&MPPSTA.VECT.TEM1')
    call jedetr('&&MPPSTA.VECT.TEM2')
    call jedetr('&&MPPSTA.VECT.TEM3')
!
!     %---------------%
!     | END OF ENDSTA |
!     %---------------%
!
end subroutine
