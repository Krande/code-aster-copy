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
subroutine mltdca(nbloc, lgbloc, ncbloc, decal, seq, &
                  nbsn, nbnd, supnd, adress, global, &
                  lgsn, factol, factou, sm, x, &
                  invp, perm, ad, trav, typsym)
! person_in_charge: olivier.boiteau at edf.fr
!     VERSION COMPLEXE DE MLTDRA
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/sspmvc.h"
#include "blas/zgemv.h"
!
    integer(kind=8) :: nbsn, nbnd, nbloc, lgbloc(nbsn), ncbloc(nbnd), decal(nbsn)
    integer(kind=4) :: global(*)
    integer(kind=8) :: seq(nbsn), supnd(nbsn+1), lgsn(nbsn)
    integer(kind=8) :: adress(nbsn+1), invp(nbnd), perm(nbnd), ad(nbnd)
    integer(kind=8) :: typsym
    complex(kind=8) :: sm(nbnd), x(nbnd), trav(nbnd)
    character(len=24) :: factol, factou, factor
    integer(kind=8) :: il, k0
    integer(kind=8) :: ib, nc, isnd, long, l, i, ndj
!
    integer(kind=8) :: seuin, seuik
    parameter(seuin=1500, seuik=300)
    integer(kind=8) :: lda, nn, kk
    integer(kind=8) :: deb1, incx, incy
    integer(kind=8) :: sni, k, j, deb, fin, adfac, ndk, gj, debndk, ifac
    complex(kind=8) :: s, alpha, beta
    character(len=1) :: tra
    blas_int :: b_incx, b_incy, b_lda, b_m, b_n
!
    call jemarq()
!
    k0 = 0
    tra = 'N'
    alpha = dcmplx(-1.d0, 0.d0)
    beta = dcmplx(1.d0, 0.d0)
    incx = 1
    incy = 1
    do j = 1, nbnd
        x(invp(j)) = sm(j)
    end do
    do j = 1, nbnd
        sm(j) = x(j)
    end do
!     DESCENTE  L * Y = B
    isnd = 0
    do ib = 1, nbloc
        call jeveuo(jexnum(factol, ib), 'L', ifac)
        adfac = ifac-1
        do nc = 1, ncbloc(ib)
            isnd = isnd+1
            sni = seq(isnd)
            long = adress(sni+1)-adress(sni)
            l = lgsn(sni)
            k = 1
            do i = adress(sni), adress(sni+1)-1
                trav(k) = x(global(i))
                k = k+1
            end do
            ad(1) = decal(sni)
            ndj = supnd(sni)-1
            do j = 1, l-1
                ndj = ndj+1
!                 RANGEMENT DU TERME DIAGONAL
                sm(ndj) = zc(ifac-1+ad(j))
!
                k = 1
                do i = j+1, l
                    trav(i) = trav(i)-zc(ifac-1+ad(j)+k)*trav(j)
                    k = k+1
                end do
!MODIF POUR SGEMV   AD(J+1) = AD(J) + LONG - J + 1
!MODIF POUR SGEMV AD(J) = AD(J) + L - J + 1
                ad(j+1) = ad(j)+long+1
                ad(j) = ad(j)+l-j+1
            end do
            ndj = ndj+1
!                 RANGEMENT DU TERME DIAGONAL
            sm(ndj) = zc(ifac-1+ad(l))
            ad(l) = ad(l)+1
            if (long .gt. l) then
!                  CALL SSPMVC(LONG-L,L,ZC(IFAC),AD,TRAV,
!     +                         TRAV(L+1))
                nn = long-l
                kk = l
                lda = long
                if (nn .lt. seuin .or. kk .lt. seuik) then
                    call sspmvc(nn, kk, zc(ifac), ad, trav, &
                                trav(l+1))
                else
                    b_lda = to_blas_int(lda)
                    b_m = to_blas_int(nn)
                    b_n = to_blas_int(kk)
                    b_incx = to_blas_int(incx)
                    b_incy = to_blas_int(incy)
                    call zgemv(tra, b_m, b_n, alpha, zc(ifac+ad(1)-1), &
                               b_lda, trav, b_incx, beta, trav(l+1), &
                               b_incy)
                end if
            end if
            k = 1
            do i = adress(sni), adress(sni+1)-1
                x(global(i)) = trav(k)
                k = k+1
            end do
        end do
        call jelibe(jexnum(factol, ib))
    end do
!
    if (typsym .ne. 0) then
        factor = factol
        deb1 = 1
    else
        factor = factou
        deb1 = nbnd
    end if
!=======================================================================
!     D * Z = Y
    do j = deb1, nbnd
        x(j) = x(j)/sm(j)
    end do
!=======================================================================
!     REMONTEE  U * X = Z
    isnd = nbsn+1
    do ib = nbloc, 1, -1
        call jeveuo(jexnum(factor, ib), 'L', ifac)
        if (ib .ne. nbloc) then
            adfac = lgbloc(ib)+ifac
        else
            adfac = lgbloc(ib)+ifac-lgsn(nbsn)
        end if
        do nc = 1, ncbloc(ib)
            isnd = isnd-1
            sni = seq(isnd)
            l = lgsn(sni)
            fin = adress(sni+1)-1
            if (sni .eq. nbsn) then
                debndk = supnd(sni+1)-2
                deb = adress(sni)+lgsn(sni)-1
                il = l-1
            else
                deb = adress(sni)+lgsn(sni)
                debndk = supnd(sni+1)-1
                il = l
            end if
            if (l .gt. 1) then
                k = 1
                do i = adress(sni), adress(sni+1)-1
                    trav(k) = x(global(i))
                    k = k+1
                end do
                k0 = k
            end if
            do ndk = debndk, supnd(sni), -1
                s = 0.d0
                if (l .gt. 1) then
                    k = k0
                    do j = fin, deb, -1
                        adfac = adfac-1
                        k = k-1
                        s = s+zc(adfac)*trav(k)
                    end do
                    deb = deb-1
                    adfac = adfac-1
                    trav(il) = trav(il)-s
                    if (typsym .eq. 0) trav(il) = trav(il)/zc(adfac)
! DECALAGE  POUR SGEMV
                    adfac = adfac-(ndk-supnd(sni))
                else
                    k = k0
                    do j = fin, deb, -1
                        gj = global(j)
                        adfac = adfac-1
                        k = k-1
                        s = s+zc(adfac)*x(gj)
                    end do
                    deb = deb-1
                    adfac = adfac-1
                    x(ndk) = x(ndk)-s
                    if (typsym .eq. 0) x(ndk) = x(ndk)/zc(adfac)
! DECALAGE  POUR SGEMV
                    adfac = adfac-(ndk-supnd(sni))
                end if
                il = il-1
            end do
            if (l .gt. 1) then
                k = 1
                do i = adress(sni), adress(sni+1)-1
                    x(global(i)) = trav(k)
                    k = k+1
                end do
            end if
        end do
        call jelibe(jexnum(factor, ib))
    end do
!     ON RANGE DANS SM  LA SOLUTION DANS LA NUMEROTATION INITIALE
    do j = 1, nbnd
        sm(perm(j)) = x(j)
    end do
    call jedema()
end subroutine
