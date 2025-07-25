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
subroutine indlia(modgen, seliai, nindep, nbddl, sst, &
                  sizlia)
!    M. CORUS     DATE 25/01/10
!-----------------------------------------------------------------------
!
!  BUT:      < CONSTRUIRE LE SOUS ESPACE POUR L'ELIMINATION >
!
!  ON IMPOSE DES LIAISONS DU TYPE C.Q=0. ON CONSTRUIT UN BASE R DU NOYAU
!  DE C. LES DDL GENERALISEES Y VERIFIENT DONC NATURELLEMENT C.R.Y=0, ET
!  Q=R.Y
!
!-----------------------------------------------------------------------
!  MODGEN  /I/ : NOM DU CONCEPT DE MODELE GENERALISE
!  SELIAI  /O/ : BASE DU NOYAU DES EQUATIONS DE LIAISON
!  NINDEP  /O/ : NOMBRE DE LIAISONS INDEPENDANTES ENTRE LES SOUS
!                  STRUCTURES
!  NBDDL   /O/ : NOMBRE DE DDL IMPLIQUES DANS LES LIAISONS
!  SST     /O/ : VECTEUR CONTENANT LES NOMS DES SOUS STRUCTURES
!  SIZLIA  /O/ : VECTEUR CONTENANT LE NB DE DDL DE CHAQUE SOUS STRUCTURE
!-----------------------------------------------------------------------
    implicit none
!
#include "jeveux.h"
#include "asterc/matfpe.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
#include "blas/dgeqp3.h"
#include "blas/dorgqr.h"
!
!-----------------------------------------------------------------------
    integer(kind=8) :: nindep, nbddl
    character(len=8) :: modgen
    character(len=24) :: seliai, sizlia, sst
!
    integer(kind=8) :: i1, j1, k1, l1, m1, nbsst, nblia, ne, nd1, nd2, lnoli1, lnoli2
    integer(kind=8) :: nedec, nd1deq, nd2deq, nbeqt, inds, lds, ldelia, llprof, lknoms
    integer(kind=8) :: lmalia, lsilia, lwork, jwork, lselia, ltau, neq, k, jjpvt, rang
    integer(kind=4) :: info
!
    character(len=8) :: int1, int2
    character(len=24) :: deflia, fprofl, nomsst, nomlia, matlia
    real(kind=8) :: eps, swork(1), x1, x2, x2prev
    blas_int :: b_lda, b_lwork, b_m, b_n
    blas_int :: b_k
    parameter(eps=2.3d-16)
!-----------------------------------------------------------------------
    call jemarq()
!
!----------------------------------------------C
!--                                          --C
!-- INITIALISATION DES DIFFERENTES GRANDEURS --C
!--                                          --C
!----------------------------------------------C
!
    deflia = modgen//'      .MODG.LIDF'
    fprofl = modgen//'      .MODG.LIPR'
    nomsst = modgen//'      .MODG.SSNO'
    nomlia = modgen//'      .MODG.LIMA'
!
!   -- NOMBRE DE SOUS STRUCTURES
    call jelira(nomsst, 'NOMMAX', nbsst)
!
!   -- NOMBRE D'INTERFACES
    call jelira(deflia, 'NMAXOC', nblia)
!
!   -- LISTE DES NOMS DES SOUS-STRUCTURES
    call wkvect(sst, 'G V K8', nbsst, lknoms)
!
!   -- LISTE DES TAILLES DES SOUS-STRUCTURES
    call wkvect(sizlia, 'G V I', nbsst, lsilia)
!
!   -- ON DETERMINE LA TAILLE DE LA MATRICE DES LIAISONS
    call jeveuo(fprofl, 'L', llprof)
    nbddl = 0
    nbeqt = 0
!
    do i1 = 1, nblia
!
!      -- NOMBRE D'EQUATIONS
        ne = zi(llprof+(i1-1)*9)
        nd1 = zi(llprof+(i1-1)*9+1)
        nd2 = zi(llprof+(i1-1)*9+4)
        nbeqt = nbeqt+ne
!
!       -- NOM DES SOUS STRUCTURES ET TAILLES
        call jeveuo(jexnum(deflia, i1), 'L', ldelia)
        int1 = zk8(ldelia)
        int2 = zk8(ldelia+2)
        if (i1 .eq. 1) then
            zk8(lknoms) = int1
            zk8(lknoms+1) = int2
            zi(lsilia) = nd1
            zi(lsilia+1) = nd2
            k1 = 2
            nbddl = nbddl+nd1+nd2
        else
            l1 = 0
            m1 = 0
            do j1 = 1, k1
                if (int1 .eq. zk8(lknoms+j1-1)) then
                    l1 = 1
                end if
                if (int2 .eq. zk8(lknoms+j1-1)) then
                    m1 = 1
                end if
            end do
            if (l1 .eq. 0) then
                zk8(lknoms+k1) = int1
                zi(lsilia+k1) = nd1
                nbddl = nbddl+nd1
                k1 = k1+1
            end if
            if (m1 .eq. 0) then
                zk8(lknoms+k1) = int2
                zi(lsilia+k1) = nd2
                nbddl = nbddl+nd2
                k1 = k1+1
            end if
        end if
    end do
!
!--------------------------------------------------------------C
!--                                                          --C
!-- ALLOCATION DE LA MATRICE L CONTENANT TOUTES LES LIAISONS --C
!--   ON CONSTRUIT SA TRANSPOSEE L^T, DIRECTEMENT UTILISABLE --C
!--   POUR LA DECOMPOSITION QR (CONSTRUCTION DU NOYAU DE L)  --C
!--                                                          --C
!--------------------------------------------------------------C
    matlia = '&&INDLIA.MATR_LIAI'
!
!-- MATRICE CARREE NEQ*NBDDL
    neq = max(nbddl, nbeqt)
    call wkvect(matlia, 'V V R', neq*nbddl, lmalia)
!
!-- ON PARCOURS LES INTERFACES POUR LA REMPLIR
!
    nedec = 0
    do k1 = 1, nblia
        ne = zi(llprof+(k1-1)*9)
        nd1 = zi(llprof+(k1-1)*9+1)
        nd2 = zi(llprof+(k1-1)*9+4)
!
!       -- RECHERCHE DE LA POSITION DE LA SOUS MATRICE DE LA
!       -- LIAISON COURANTE DANS LA MATRICE GLOBALE
!
        nd1deq = 0
        nd2deq = 0
        call jeveuo(jexnum(deflia, k1), 'L', ldelia)
        int1 = zk8(ldelia)
        int2 = zk8(ldelia+2)
        do i1 = 1, nbsst
            if (int1 .eq. zk8(lknoms+i1-1)) then
                do j1 = 1, i1-1
                    nd1deq = nd1deq+zi(lsilia+j1-1)
                end do
            end if
            if (int2 .eq. zk8(lknoms+i1-1)) then
                do j1 = 1, i1-1
                    nd2deq = nd2deq+zi(lsilia+j1-1)
                end do
            end if
        end do
!
!       -- REMPLISSAGE DE LA SOUS MATRICE POUR LA LIAISON K1
!
        call jeveuo(jexnum(nomlia, (k1-1)*3+1), 'L', lnoli1)
        call jeveuo(jexnum(nomlia, (k1-1)*3+2), 'L', lnoli2)
!
        do j1 = 1, ne
            do i1 = 1, nd1
                if (abs(zr(lnoli1+(i1-1)*ne+j1-1)) .gt. eps) then
                    zr(lmalia+(j1-1+nedec)*nbddl+i1-1+nd1deq) = zr(lnoli1+(i1-1)*ne+j1-1)
                else
                    zr(lmalia+(j1-1+nedec)*nbddl+i1-1+nd1deq) = 0.d0
                end if
            end do
!
            do i1 = 1, nd2
                if (abs(zr(lnoli2+(i1-1)*ne+j1-1)) .gt. eps) then
                    zr(lmalia+(j1-1+nedec)*nbddl+i1-1+nd2deq) = zr(lnoli2+(i1-1)*ne+j1-1)
                else
                    zr(lmalia+(j1-1+nedec)*nbddl+i1-1+nd2deq) = 0.d0
                end if
            end do
        end do
!
        nedec = nedec+ne
    end do
!
!
!-------------------------------------------------------------C
!--                                                         --C
!-- QR DE LA MATRICE POUR DETERMINER LA TAILLE DE SON NOYAU --C
!--                                                         --C
!-------------------------------------------------------------C
!
    lds = int(min(neq, nbddl))
!
    call wkvect('&&INDLIA.TAU', 'V V R', neq, ltau)
!
!   -- desactivation du test fpe
    call matfpe(-1)
!
!
!   -- Decomposition QR  avec dgeqp3 (avec pivotage) :
!   ----------------------------------------------------
    call wkvect('&&INDLIA.JPVT', 'V V S', neq, jjpvt)
    lwork = -1
    b_lda = to_blas_int(nbddl)
    b_m = to_blas_int(nbddl)
    b_n = to_blas_int(neq)
    b_lwork = to_blas_int(lwork)
    call dgeqp3(b_m, b_n, zr(lmalia), b_lda, zi4(jjpvt), &
                zr(ltau), swork, b_lwork, info)
    ASSERT(info .eq. 0)
    lwork = int(swork(1))
    call wkvect('&&MATR_QR_WORK', 'V V R', lwork, jwork)
!
    do k = 1, nbddl*neq
    end do
    b_lda = to_blas_int(nbddl)
    b_m = to_blas_int(nbddl)
    b_n = to_blas_int(neq)
    b_lwork = to_blas_int(lwork)
    call dgeqp3(b_m, b_n, zr(lmalia), b_lda, zi4(jjpvt), &
                zr(ltau), zr(jwork), b_lwork, info)
    ASSERT(info .eq. 0)
!
!
!   -- Grace a dgeqp3, on sait que la diagonale de R est "decroissante" :
!      on cherche a determiner le rang de matlia on regardant la diagonale de R
!   ----------------------------------------------------------------------------
    x1 = abs(zr(lmalia))
    x2prev = x1
    do i1 = 2, neq
        x2 = abs(zr(lmalia-1+(i1-1)*nbddl+i1))
!
!       -- soit x2 est tres petit devant x1 :
        if ((x2/x1) .lt. 1.d-12) then
            rang = i1-1
            exit
        end if
!
!       -- soit x2 est beaucoup plus petit que le terme precedent :
        if ((x2/x2prev) .lt. 1.d-5) then
            rang = i1-1
            exit
        end if
        x2prev = x2
    end do
    nindep = nbddl-rang
!
!
    write (6, *) '--------'
    write (6, *) ' '
    write (6, *) '+++', nbddl, ' DEGRES DE LIBERTE AU TOTAL'
    write (6, *) '+++', nbeqt, ' CONTRAINTES CINEMATIQUES.'
    write (6, *) ' '
    write (6, *) 'ON A TROUVE', nindep, ' RELATIONS INDEPENDANTES.'
    write (6, *) ' '
    write (6, *) '--------'
!
!
!   -- Il ne peut pas y avoir moins de ddl independants que le nombre
!   --  total de ddl moins le nombre de contraintes.
!   ----------------------------------------------------------------------
    ASSERT(nindep .ge. (nbddl-nbeqt))
!
!
!
!   -- Reconstruction de la matrice q  (dorgqr):
!   ---------------------------------------------
    b_lda = to_blas_int(nbddl)
    b_m = to_blas_int(nbddl)
    b_n = to_blas_int(neq)
    b_k = to_blas_int(nbddl)
    b_lwork = to_blas_int(-1)
    call dorgqr(b_m, b_n, b_k, zr(lmalia), b_lda, &
                zr(ltau), swork(1), b_lwork, info)
    ASSERT(info .eq. 0)
    if (swork(1) .gt. lwork) then
        lwork = int(swork(1))
        call jedetr('&&MATR_QR_WORK')
        call wkvect('&&MATR_QR_WORK', 'V V R', lwork, jwork)
    end if
    b_lda = to_blas_int(nbddl)
    b_m = to_blas_int(nbddl)
    b_n = to_blas_int(neq)
    b_k = to_blas_int(nbddl)
    b_lwork = to_blas_int(lwork)
    call dorgqr(b_m, b_n, b_k, zr(lmalia), b_lda, &
                zr(ltau), zr(jwork), b_lwork, info)
    ASSERT(info .eq. 0)
!
!   -- reactivation du test fpe
    call matfpe(1)
!
!
!   -- Construction du sous espace :
!   ---------------------------------
    call wkvect(seliai, 'G V R', nindep*nbddl, lselia)
    do j1 = 1, nindep
        inds = nbddl+1-j1
        do i1 = 1, nbddl
            zr(lselia+(j1-1)*nbddl+i1-1) = zr(lmalia+(inds-1)*nbddl+i1-1)
        end do
    end do
!
!
!   -- destruction des objets de travail :
!   -------------------------------------------
    call jedetr(matlia)
    call jedetr('&&MATR_QR_WORK')
    call jedetr('&&INDLIA.TAU')
    call jedetr('&&INDLIA.JPVT')
!
    call jedema()
!
end subroutine
