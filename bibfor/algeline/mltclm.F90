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
subroutine mltclm(nb, n, p, front, adper, &
                  t1, ad, eps, ier, c)
! person_in_charge: olivier.boiteau at edf.fr
    implicit none
#include "asterfort/mltcld.h"
#include "asterfort/mltclj.h"
#include "blas/zgemv.h"
    integer(kind=8) :: n, p, adper(*), ad(*), ier, nb
    real(kind=8) :: eps
    complex(kind=8) :: front(*), t1(*), c(nb, nb, *), alpha, beta
    integer(kind=8) :: i, kb, adk, adki, decal, l
    integer(kind=8) :: m, ll, k, ind, ia, j, restp, npb
    integer(kind=8) :: incx, incy
    character(len=1) :: tra
    blas_int :: b_incx, b_incy, b_lda, b_m, b_n
    npb = p/nb
    restp = p-(nb*npb)
    ll = n
    tra = 'N'
    alpha = dcmplx(-1.d0, 0.d0)
    beta = dcmplx(1.d0, 0.d0)
    incx = 1
    incy = 1
!
    do kb = 1, npb
!     K : INDICE (DANS LA MATRICE FRONTALE ( DE 1 A P)),
!     DE LA PREMIERE COLONNE DU BLOC
        k = nb*(kb-1)+1
        adk = adper(k)
!     BLOC DIAGONAL
        call mltcld(nb, front(adk), adper, t1, ad, &
                    eps, ier)
        if (ier .gt. 0) goto 999
!
!     NORMALISATION DES BLOCS SOUS LE BLOC DIAGONAL
!
        ll = ll-nb
        ia = adk+nb
        do i = 1, nb
            ind = ia+n*(i-1)
            if (i .gt. 1) then
                do l = 1, i-1
                    t1(l) = front(adper(k+l-1))*front(n*(k+l-2)+k+i-1)
                end do
            end if
            b_lda = to_blas_int(n)
            b_m = to_blas_int(ll)
            b_n = to_blas_int(i-1)
            b_incx = to_blas_int(incx)
            b_incy = to_blas_int(incy)
            call zgemv(tra, b_m, b_n, alpha, front(ia), &
                       b_lda, t1, b_incx, beta, front(ind), &
                       b_incy)
            adki = adper(k+i-1)
            do j = 1, ll
                front(ind) = front(ind)/front(adki)
                ind = ind+1
            end do
        end do
!
        decal = kb*nb
        ll = n-decal
        m = p-decal
        ind = adper(k+nb)
        call mltclj(nb, n, ll, m, k, &
                    decal, front, front(ind), adper, t1, &
                    c)
    end do
!     COLONNES RESTANTES
    if (restp .gt. 0) then
!     K : INDICE (DANS LA MATRICE FRONTALE ( DE 1 A P)),
!     DE LA PREMIERE COLONNE DU BLOC
        kb = npb+1
        k = nb*npb+1
        adk = adper(k)
!     BLOC DIAGONAL
        call mltcld(restp, front(adk), adper, t1, ad, &
                    eps, ier)
        if (ier .gt. 0) goto 999
!
!     NORMALISATION DES BLOCS SOUS LE BLOC DIAGONAL
!
        ll = n-p
        ia = adk+restp
        do i = 1, restp
            ind = ia+n*(i-1)
            if (i .gt. 1) then
                do l = 1, i-1
                    t1(l) = front(adper(k+l-1))*front(n*(k+l-2)+k+i-1)
                end do
            end if
            b_lda = to_blas_int(n)
            b_m = to_blas_int(ll)
            b_n = to_blas_int(i-1)
            b_incx = to_blas_int(incx)
            b_incy = to_blas_int(incy)
            call zgemv(tra, b_m, b_n, alpha, front(ia), &
                       b_lda, t1, b_incx, beta, front(ind), &
                       b_incy)
            adki = adper(k+i-1)
            do j = 1, ll
                front(ind) = front(ind)/front(adki)
                ind = ind+1
            end do
        end do
!
    end if
999 continue
    if (ier .gt. 0) ier = ier+nb*(kb-1)
end subroutine
