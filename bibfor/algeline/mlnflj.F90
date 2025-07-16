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
subroutine mlnflj(nb, n, ll, m, it, &
                  p, frontl, frontu, frnl, frnu, &
                  adper, travl, travu, cl, cu)
! person_in_charge: olivier.boiteau at edf.fr
    use superv_module
    implicit none
! aslint: disable=C1513
#include "blas/dgemm.h"
    integer(kind=8) :: n, p, adper(*)
    real(kind=8) :: frontl(*), frnl(*), frontu(*), frnu(*)
    integer(kind=8) :: nb, decal, add, ind, nmb, i, j, kb, ia, ib, nlb, ll
    character(len=1) :: tra, trb
    integer(kind=8) :: m, k, i1, it, j1, restm, restl, nbl
    integer(kind=8) :: nproc, numpro
    real(kind=8) :: travl(p, nb, *), travu(p, nb, *)
    real(kind=8) :: cl(nb, nb, *), cu(nb, nb, *), alpha, beta
    blas_int :: b_k, b_lda, b_ldb, b_ldc, b_m, b_n
    tra = 'N'
    trb = 'N'
    alpha = -1.d0
    beta = 0.d0
    nbl = p-it+1
    nmb = m/nb
    nlb = ll/nb
    restm = m-(nb*nmb)
    restl = ll-(nb*nlb)
    decal = adper(p+1)-1
    nproc = asthread_getmax()
    numpro = 1
    if (nmb .ge. nproc) then
        !$OMP PARALLEL DO DEFAULT(PRIVATE) &
        !$OMP SHARED(N,M,P,NMB,NBL,NLB,NB,RESTM,RESTL,CL,CU,TRA,TRB,ALPHA,BETA) &
        !$OMP SHARED(FRONTL,FRONTU,ADPER,DECAL,FRNL,FRNU,TRAVL,TRAVU,IT) &
        !$OMP SCHEDULE(STATIC,1)
        do kb = 1, nmb
            numpro = asthread_getnum()+1
!     K : INDICE DE COLONNE DANS LA MATRICE FRONTALE (ABSOLU DE 1 A N)
            k = nb*(kb-1)+1+p
            do i = it, p
                add = n*(i-1)+k
                do j = 1, nb
                    travl(i, j, numpro) = frontl(add)
                    travu(i, j, numpro) = frontu(add)
                    add = add+1
                end do
            end do
!     BLOC DIAGONAL
!
!     SOUS LE BLOC DIAGONAL
!     2EME ESSAI : DES PRODUITS DE LONGUEUR NB
!
!
            do ib = kb, nlb
                ia = n*(it-1)+k+nb*(ib-kb)
                b_ldc = to_blas_int(nb)
                b_ldb = to_blas_int(p)
                b_lda = to_blas_int(n)
                b_m = to_blas_int(nb)
                b_n = to_blas_int(nb)
                b_k = to_blas_int(nbl)
                call dgemm(tra, trb, b_m, b_n, b_k, &
                           alpha, frontl(ia), b_lda, travu(it, 1, numpro), b_ldb, &
                           beta, cl(1, 1, numpro), b_ldc)
                b_ldc = to_blas_int(nb)
                b_ldb = to_blas_int(p)
                b_lda = to_blas_int(n)
                b_m = to_blas_int(nb)
                b_n = to_blas_int(nb)
                b_k = to_blas_int(nbl)
                call dgemm(tra, trb, b_m, b_n, b_k, &
                           alpha, frontu(ia), b_lda, travl(it, 1, numpro), b_ldb, &
                           beta, cu(1, 1, numpro), b_ldc)
!     RECOPIE
!
!
                do i = 1, nb
                    i1 = i-1
!              IND = ADPER(K +I1) - DECAL  + NB*(IB-KB-1) +NB - I1
                    if (ib .eq. kb) then
                        j1 = i
                        ind = adper(k+i1)-decal
                    else
                        j1 = 1
                        ind = adper(k+i1)-decal+nb*(ib-kb)-i1
                    end if
                    do j = j1, nb
                        frnl(ind) = frnl(ind)+cl(j, i, numpro)
                        frnu(ind) = frnu(ind)+cu(j, i, numpro)
                        ind = ind+1
                    end do
                end do
            end do
            if (restl .gt. 0) then
                ib = nlb+1
                ia = n*(it-1)+k+nb*(ib-kb)
                b_ldc = to_blas_int(nb)
                b_ldb = to_blas_int(p)
                b_lda = to_blas_int(n)
                b_m = to_blas_int(restl)
                b_n = to_blas_int(nb)
                b_k = to_blas_int(nbl)
                call dgemm(tra, trb, b_m, b_n, b_k, &
                           alpha, frontl(ia), b_lda, travu(it, 1, numpro), b_ldb, &
                           beta, cl(1, 1, numpro), b_ldc)
                b_ldc = to_blas_int(nb)
                b_ldb = to_blas_int(p)
                b_lda = to_blas_int(n)
                b_m = to_blas_int(restl)
                b_n = to_blas_int(nb)
                b_k = to_blas_int(nbl)
                call dgemm(tra, trb, b_m, b_n, b_k, &
                           alpha, frontu(ia), b_lda, travl(it, 1, numpro), b_ldb, &
                           beta, cu(1, 1, numpro), b_ldc)
!           RECOPIE
!
!
                do i = 1, nb
                    i1 = i-1
                    j1 = 1
                    ind = adper(k+i1)-decal+nb*(ib-kb)-i1
                    do j = j1, restl
                        frnl(ind) = frnl(ind)+cl(j, i, numpro)
                        frnu(ind) = frnu(ind)+cu(j, i, numpro)
                        ind = ind+1
                    end do
                end do
            end if
        end do
    else
        do kb = 1, nmb
!     K : INDICE DE COLONNE DANS LA MATRICE FRONTALE (ABSOLU DE 1 A N)
            k = nb*(kb-1)+1+p
            do i = it, p
                add = n*(i-1)+k
                do j = 1, nb
                    travl(i, j, 1) = frontl(add)
                    travu(i, j, 1) = frontu(add)
                    add = add+1
                end do
            end do
!     BLOC DIAGONAL
!
!     SOUS LE BLOC DIAGONAL
!     2EME ESSAI : DES PRODUITS DE LONGUEUR NB
!
!
            do ib = kb, nlb
                ia = n*(it-1)+k+nb*(ib-kb)
                b_ldc = to_blas_int(nb)
                b_ldb = to_blas_int(p)
                b_lda = to_blas_int(n)
                b_m = to_blas_int(nb)
                b_n = to_blas_int(nb)
                b_k = to_blas_int(nbl)
                call dgemm(tra, trb, b_m, b_n, b_k, &
                           alpha, frontl(ia), b_lda, travu(it, 1, 1), b_ldb, &
                           beta, cl(1, 1, 1), b_ldc)
                b_ldc = to_blas_int(nb)
                b_ldb = to_blas_int(p)
                b_lda = to_blas_int(n)
                b_m = to_blas_int(nb)
                b_n = to_blas_int(nb)
                b_k = to_blas_int(nbl)
                call dgemm(tra, trb, b_m, b_n, b_k, &
                           alpha, frontu(ia), b_lda, travl(it, 1, 1), b_ldb, &
                           beta, cu(1, 1, 1), b_ldc)
!     RECOPIE
!
!
                do i = 1, nb
                    i1 = i-1
                    if (ib .eq. kb) then
                        j1 = i
                        ind = adper(k+i1)-decal
                    else
                        j1 = 1
                        ind = adper(k+i1)-decal+nb*(ib-kb)-i1
                    end if
                    do j = j1, nb
                        frnl(ind) = frnl(ind)+cl(j, i, 1)
                        frnu(ind) = frnu(ind)+cu(j, i, 1)
                        ind = ind+1
                    end do
                end do
            end do
            if (restl .gt. 0) then
                ib = nlb+1
                ia = n*(it-1)+k+nb*(ib-kb)
                b_ldc = to_blas_int(nb)
                b_ldb = to_blas_int(p)
                b_lda = to_blas_int(n)
                b_m = to_blas_int(restl)
                b_n = to_blas_int(nb)
                b_k = to_blas_int(nbl)
                call dgemm(tra, trb, b_m, b_n, b_k, &
                           alpha, frontl(ia), b_lda, travu(it, 1, 1), b_ldb, &
                           beta, cl(1, 1, 1), b_ldc)
                b_ldc = to_blas_int(nb)
                b_ldb = to_blas_int(p)
                b_lda = to_blas_int(n)
                b_m = to_blas_int(restl)
                b_n = to_blas_int(nb)
                b_k = to_blas_int(nbl)
                call dgemm(tra, trb, b_m, b_n, b_k, &
                           alpha, frontu(ia), b_lda, travl(it, 1, 1), b_ldb, &
                           beta, cu(1, 1, 1), b_ldc)
!           RECOPIE
!
!
                do i = 1, nb
                    i1 = i-1
!              IND = ADPER(K +I1) - DECAL  + NB*(IB-KB-1) +NB - I1
                    j1 = 1
                    ind = adper(k+i1)-decal+nb*(ib-kb)-i1
                    do j = j1, restl
                        frnl(ind) = frnl(ind)+cl(j, i, 1)
                        frnu(ind) = frnu(ind)+cu(j, i, 1)
                        ind = ind+1
                    end do
                end do
            end if
        end do
    end if
    if (restm .gt. 0) then
        kb = 1+nmb
!     K : INDICE DE COLONNE DANS LA MATRICE FRONTALE (ABSOLU DE 1 A N)
        k = nb*(kb-1)+1+p
        do i = it, p
            add = n*(i-1)+k
            do j = 1, restm
                travl(i, j, 1) = frontl(add)
                travu(i, j, 1) = frontu(add)
                add = add+1
            end do
        end do
!     BLOC DIAGONAL
!
!     SOUS LE BLOC DIAGONAL
!     2EME ESSAI : DES PRODUITS DE LONGUEUR NB
!
        do ib = kb, nlb
            ia = n*(it-1)+k+nb*(ib-kb)
            b_ldc = to_blas_int(nb)
            b_ldb = to_blas_int(p)
            b_lda = to_blas_int(n)
            b_m = to_blas_int(nb)
            b_n = to_blas_int(restm)
            b_k = to_blas_int(nbl)
            call dgemm(tra, trb, b_m, b_n, b_k, &
                       alpha, frontl(ia), b_lda, travu(it, 1, 1), b_ldb, &
                       beta, cl(1, 1, 1), b_ldc)
            b_ldc = to_blas_int(nb)
            b_ldb = to_blas_int(p)
            b_lda = to_blas_int(n)
            b_m = to_blas_int(nb)
            b_n = to_blas_int(restm)
            b_k = to_blas_int(nbl)
            call dgemm(tra, trb, b_m, b_n, b_k, &
                       alpha, frontu(ia), b_lda, travl(it, 1, 1), b_ldb, &
                       beta, cu(1, 1, 1), b_ldc)
!     RECOPIE
!
!
            do i = 1, restm
                i1 = i-1
!     IND = ADPER(K +I1) - DECAL  + NB*(IB-KB-1) +NB - I1
                if (ib .eq. kb) then
                    j1 = i
                    ind = adper(k+i1)-decal
                else
                    j1 = 1
                    ind = adper(k+i1)-decal+nb*(ib-kb)-i1
                end if
                do j = j1, nb
                    frnl(ind) = frnl(ind)+cl(j, i, 1)
                    frnu(ind) = frnu(ind)+cu(j, i, 1)
                    ind = ind+1
                end do
            end do
        end do
        if (restl .gt. 0) then
            ib = nlb+1
            ia = n*(it-1)+k+nb*(ib-kb)
            b_ldc = to_blas_int(nb)
            b_ldb = to_blas_int(p)
            b_lda = to_blas_int(n)
            b_m = to_blas_int(restl)
            b_n = to_blas_int(restm)
            b_k = to_blas_int(nbl)
            call dgemm(tra, trb, b_m, b_n, b_k, &
                       alpha, frontl(ia), b_lda, travu(it, 1, 1), b_ldb, &
                       beta, cl(1, 1, 1), b_ldc)
            b_ldc = to_blas_int(nb)
            b_ldb = to_blas_int(p)
            b_lda = to_blas_int(n)
            b_m = to_blas_int(restl)
            b_n = to_blas_int(restm)
            b_k = to_blas_int(nbl)
            call dgemm(tra, trb, b_m, b_n, b_k, &
                       alpha, frontu(ia), b_lda, travl(it, 1, 1), b_ldb, &
                       beta, cu(1, 1, 1), b_ldc)
!     RECOPIE
!
!
            do i = 1, restm
                i1 = i-1
!     IND = ADPER(K +I1) - DECAL  + NB*(IB-KB-1) +NB - I1
                j1 = 1
                ind = adper(k+i1)-decal+nb*(ib-kb)-i1
                do j = j1, restl
                    frnl(ind) = frnl(ind)+cl(j, i, 1)
                    frnu(ind) = frnu(ind)+cu(j, i, 1)
                    ind = ind+1
                end do
            end do
        end if
    end if
end subroutine
