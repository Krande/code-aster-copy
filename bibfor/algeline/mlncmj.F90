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
subroutine mlncmj(nb, n, p, frontl, frontu, &
                  frnl, frnu, adper, t1, t2, &
                  cl, cu)
! person_in_charge: olivier.boiteau at edf.fr
!
!     VERSION AVEC APPEL A DGEMV POUR LES PRODUITS MATRICE-VECTEUR
!     AU DELA D' UN CERTAIN SEUIL
!     DGEMV EST APPEL A T1ERS LA FONCTION C DGEMW POUR CAR DGEMV
!     NECESSITE  DES ARGUMENTS ENTIER INTEGER*4 REFUSES PAR ASTER
!
    use superv_module
    implicit none
! aslint: disable=C1513
#include "blas/zgemm.h"
    integer(kind=8) :: n, p, adper(*), restm, decal
    complex(kind=8) :: frontl(*), frontu(*), frnl(*), frnu(*)
    integer(kind=8) :: nmb
    character(len=1) :: tra, trb
    integer(kind=8) :: i1, j1, k, m, it, nb, numprc
    complex(kind=8) :: t1(p, nb, *), t2(p, nb, *), alpha, beta
    complex(kind=8) :: cl(nb, nb, *), cu(nb, nb, *)
    integer(kind=8) :: i, kb, j, ib, ia, ind, add
    blas_int :: b_k, b_lda, b_ldb, b_ldc, b_m, b_n
    m = n-p
    nmb = m/nb
    restm = m-(nb*nmb)
    decal = adper(p+1)-1
    tra = 'N'
    trb = 'N'
    alpha = dcmplx(-1.d0, 0.d0)
    beta = dcmplx(0.d0, 0.d0)
!
    !$OMP PARALLEL DO DEFAULT(PRIVATE) &
    !$OMP SHARED(N,M,P,NMB,NB,RESTM,FRONTL,FRONTU,ADPER,DECAL,FRNL,FRNU) &
    !$OMP SHARED(T1,T2,CL,CU,TRA,TRB,ALPHA,BETA) &
    !$OMP SCHEDULE(STATIC,1)
    do kb = 1, nmb
        numprc = asthread_getnum()+1
!     K : INDICE DE COLONNE DANS LA MATRICE FRONTALE (ABSOLU DE 1 A N)
        k = nb*(kb-1)+1+p
        do i = 1, p
            add = n*(i-1)+k
            do j = 1, nb
                t1(i, j, numprc) = frontu(add)
                t2(i, j, numprc) = frontl(add)
                add = add+1
            end do
        end do
!     BLOC DIAGONAL
!
!     SOUS LE BLOC DIAGONAL
!     2EME ESSAI : DES PRODUITS DE LONGUEUR NB
!
        do ib = kb, nmb
            ia = k+nb*(ib-kb)
            it = 1
            b_ldc = to_blas_int(nb)
            b_ldb = to_blas_int(p)
            b_lda = to_blas_int(n)
            b_m = to_blas_int(nb)
            b_n = to_blas_int(nb)
            b_k = to_blas_int(p)
            call zgemm(tra, trb, b_m, b_n, b_k, &
                       alpha, frontl(ia), b_lda, t1(it, 1, numprc), b_ldb, &
                       beta, cl(1, 1, numprc), b_ldc)
            b_ldc = to_blas_int(nb)
            b_ldb = to_blas_int(p)
            b_lda = to_blas_int(n)
            b_m = to_blas_int(nb)
            b_n = to_blas_int(nb)
            b_k = to_blas_int(p)
            call zgemm(tra, trb, b_m, b_n, b_k, &
                       alpha, frontu(ia), b_lda, t2(it, 1, numprc), b_ldb, &
                       beta, cu(1, 1, numprc), b_ldc)
!     RECOPIE
!
!
            do i = 1, nb
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
                    frnl(ind) = frnl(ind)+cl(j, i, numprc)
                    frnu(ind) = frnu(ind)+cu(j, i, numprc)
                    ind = ind+1
                end do
            end do
        end do
!
        if (restm .gt. 0) then
            ib = nmb+1
            ia = k+nb*(ib-kb)
            it = 1
            b_ldc = to_blas_int(nb)
            b_ldb = to_blas_int(p)
            b_lda = to_blas_int(n)
            b_m = to_blas_int(restm)
            b_n = to_blas_int(nb)
            b_k = to_blas_int(p)
            call zgemm(tra, trb, b_m, b_n, b_k, &
                       alpha, frontl(ia), b_lda, t1(it, 1, numprc), b_ldb, &
                       beta, cl(1, 1, numprc), b_ldc)
            b_ldc = to_blas_int(nb)
            b_ldb = to_blas_int(p)
            b_lda = to_blas_int(n)
            b_m = to_blas_int(restm)
            b_n = to_blas_int(nb)
            b_k = to_blas_int(p)
            call zgemm(tra, trb, b_m, b_n, b_k, &
                       alpha, frontu(ia), b_lda, t2(it, 1, numprc), b_ldb, &
                       beta, cu(1, 1, numprc), b_ldc)
!     RECOPIE
!
!
            do i = 1, nb
                i1 = i-1
!     IND = ADPER(K +I1) - DECAL  + NB*(IB-KB-1) +NB - I1
                j1 = 1
                ind = adper(k+i1)-decal+nb*(ib-kb)-i1
                do j = j1, restm
                    frnl(ind) = frnl(ind)+cl(j, i, numprc)
                    frnu(ind) = frnu(ind)+cu(j, i, numprc)
                    ind = ind+1
                end do
            end do
!
!
        end if
    end do
    !$OMP END PARALLEL DO
    numprc = 1
    if (restm .gt. 0) then
        kb = 1+nmb
!     K : INDICE DE COLONNE DANS LA MATRICE FRONTLALE (ABSOLU DE 1 A N)
        k = nb*(kb-1)+1+p
        do i = 1, p
            add = n*(i-1)+k
            do j = 1, restm
                t1(i, j, 1) = frontu(add)
                t2(i, j, 1) = frontl(add)
                add = add+1
            end do
        end do
!     BLOC DIAGONAL
!
        ib = kb
        ia = k+nb*(ib-kb)
        it = 1
        b_ldc = to_blas_int(nb)
        b_ldb = to_blas_int(p)
        b_lda = to_blas_int(n)
        b_m = to_blas_int(restm)
        b_n = to_blas_int(restm)
        b_k = to_blas_int(p)
        call zgemm(tra, trb, b_m, b_n, b_k, &
                   alpha, frontl(ia), b_lda, t1(it, 1, 1), b_ldb, &
                   beta, cl(1, 1, numprc), b_ldc)
        b_ldc = to_blas_int(nb)
        b_ldb = to_blas_int(p)
        b_lda = to_blas_int(n)
        b_m = to_blas_int(restm)
        b_n = to_blas_int(restm)
        b_k = to_blas_int(p)
        call zgemm(tra, trb, b_m, b_n, b_k, &
                   alpha, frontu(ia), b_lda, t2(it, 1, 1), b_ldb, &
                   beta, cu(1, 1, numprc), b_ldc)
!     RECOPIE
!
!
        do i = 1, restm
            i1 = i-1
!     IND = ADPER(K +I1) - DECAL  + NB*(IB-KB-1) +NB - I1
            j1 = i
            ind = adper(k+i1)-decal
!
            do j = j1, restm
!
                frnl(ind) = frnl(ind)+cl(j, i, numprc)
                frnu(ind) = frnu(ind)+cu(j, i, numprc)
                ind = ind+1
            end do
        end do
!
    end if
end subroutine
