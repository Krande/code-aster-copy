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
subroutine projmc(matras, nomres, basemo, nugene, nu, &
                  neq, nbmo)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/copmod.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mcmult.h"
#include "asterfort/mtdscr.h"
#include "asterfort/rsexch.h"
#include "asterfort/ualfva.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/zerlag.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    integer(kind=8) :: neq, nbmo
    character(len=8) :: matras, nomres, basemo
    character(len=14) :: nu
    character(len=19) :: nomsto
    character(len=14) :: nugene
!
!     CALCUL PROJECTION MATRICE COMPLEXE SUR BASE DE RITZ
!
!-----------------------------------------------------------------------
!
    integer(kind=8) :: iddeeq, nueq, ntbloc, nbloc, iaconl, jrefa, iadesc, nbj
    integer(kind=8) :: i, j, k, imatra, iblo, ldblo, n1bloc, n2bloc, hc, ldblo1, ldblo2
    integer(kind=8) :: iret
    character(len=1) :: typbase
    character(len=19) :: matr, resu
    character(len=24) :: nomcha
    real(kind=8) :: zero
    complex(kind=8) :: cbid, pij
    aster_logical :: lsym
    real(kind=8), pointer :: vbasemo(:) => null()
    complex(kind=8), pointer :: vectass2(:) => null()
    complex(kind=8), pointer :: vectass3(:) => null()
    integer(kind=8), pointer :: smde(:) => null()
    integer(kind=4), pointer :: smhc(:) => null()
    integer(kind=8), pointer :: smdi(:) => null()
    cbid = dcmplx(0.d0, 0.d0)
!-----------------------------------------------------------------------
!
    call jemarq()
!
    zero = 0.d0
    resu = ' '
    resu = nomres
    matr = matras
!
!   TYPE DE LA BASE MODALE
    call rsexch('F', basemo, 'DEPL', 1, nomcha, &
                iret)
    call jelira(nomcha(1:19)//'.VALE', 'TYPE', cval=typbase)
    if (typbase .eq. 'C') call utmess('A', 'DEFIBASEMODALE1_2')
!
    call jeveuo(nu//'.NUME.DEEQ', 'L', iddeeq)
!
    nomsto = nugene//'.SMOS'
    call jeveuo(nomsto//'.SMDE', 'L', vi=smde)
    nueq = smde(1)
    ntbloc = smde(2)
    nbloc = smde(3)
!
    call jeveuo(matr//'.REFA', 'L', jrefa)
    lsym = zk24(jrefa-1+9) .eq. 'MS'
    if (lsym) then
        call jecrec(resu//'.UALF', 'G V C', 'NU', 'DISPERSE', 'CONSTANT', &
                    nbloc)
    else
        call jecrec(resu//'.UALF', 'G V C', 'NU', 'DISPERSE', 'CONSTANT', &
                    2*nbloc)
    end if
!
    call jeecra(resu//'.UALF', 'LONMAX', ntbloc)
!
    call wkvect(resu//'.CONL', 'G V C', nueq, iaconl)
    do i = 1, nueq
        zc(iaconl+i-1) = dcmplx(1.0d0, 0.0d0)
    end do
!
    call wkvect(resu//'.REFA', 'G V K24', 20, jrefa)
    zk24(jrefa-1+11) = 'MPI_COMPLET'
    zk24(jrefa-1+1) = basemo
    zk24(jrefa-1+2) = nugene
    if (lsym) then
        zk24(jrefa-1+9) = 'MS'
    else
        zk24(jrefa-1+9) = 'MR'
    end if
    zk24(jrefa-1+10) = 'GENE'
!
    call wkvect(resu//'.DESC', 'G V I', 3, iadesc)
    zi(iadesc) = 2
    zi(iadesc+1) = nueq
!   ON TESTE LA HAUTEUR MAXIMALE DES COLONNES DE LA MATRICE
!   SI CETTE HAUTEUR VAUT 1, ON SUPPOSE QUE LE STOCKAGE EST DIAGONAL
    if (smde(2) .eq. smde(1)) then
        zi(iadesc+2) = 1
    else
        zi(iadesc+2) = 2
    end if
!
    AS_ALLOCATE(vc=vectass2, size=neq)
    AS_ALLOCATE(vc=vectass3, size=neq)
    AS_ALLOCATE(vr=vbasemo, size=nbmo*neq)
! ----- CONVERSION DE BASEMO A LA NUMEROTATION NU
    call copmod(basemo, numer=nu, bmodr=vbasemo)
!
    call mtdscr(matras)
    call jeveuo(matras//'           .&INT', 'E', imatra)
!
! --- RECUPERATION DE LA STRUCTURE DE LA MATR_ASSE_GENE
!
    call jeveuo(nomsto//'.SMDI', 'L', vi=smdi)
    call jeveuo(nomsto//'.SMHC', 'L', vi4=smhc)
!
!     -- CAS DES MATRICES SYMETRIQUES :
!     ----------------------------------
    if (lsym) then
        do iblo = 1, nbloc
!
            call jecroc(jexnum(resu//'.UALF', iblo))
            call jeveuo(jexnum(resu//'.UALF', iblo), 'E', ldblo)
!
! ------ PROJECTION DE LA MATRICE ASSEMBLEE
!
!        BOUCLE SUR LES COLONNES DE LA MATRICE ASSEMBLEE
!
            n1bloc = 1
            n2bloc = smde(1)
!
            do i = n1bloc, n2bloc
!
! --------- MISE A ZERO PARTIE IMAGINAIRE DU MODE I
!
                do k = 1, neq
                    vectass2(k) = dcmplx(vbasemo(1+(i-1)*neq+k-1), zero)
                end do
!
! --------- CALCUL PRODUIT MATRICE*MODE I
!
                call mcmult('ZERO', imatra, vectass2, vectass3, 1, &
                            .true._1)
                call zerlag(neq, zi(iddeeq), vectz=vectass3)
!
! --------- BOUCLE SUR LES INDICES VALIDES DE LA COLONNE I
!
                hc = smdi(i)
                if (i .gt. 1) hc = hc-smdi(i-1)
                do j = (i-hc+1), i
!
! ----------- PRODUIT SCALAIRE VECTASS * MODE
!
                    pij = dcmplx(zero, zero)
                    do k = 1, neq
                        pij = pij+vectass3(k)*dcmplx(vbasemo(1+(j-1)*neq+k-1), zero)
                    end do
!
! ----------- STOCKAGE DANS LE .UALF A LA BONNE PLACE (1 BLOC)
!
                    zc(ldblo+smdi(i)+j-i-1) = pij
!
                end do
            end do
            call jelibe(jexnum(resu//'.UALF', iblo))
        end do
    else
!     -- CAS DES MATRICES NON-SYMETRIQUES :
!     --------------------------------------
!
        ASSERT(nbloc .eq. 1)
        call jecroc(jexnum(resu//'.UALF', 1))
        call jecroc(jexnum(resu//'.UALF', 2))
        call jeveuo(jexnum(resu//'.UALF', 1), 'E', ldblo1)
        call jeveuo(jexnum(resu//'.UALF', 2), 'E', ldblo2)
        n1bloc = 1
        n2bloc = smde(1)
        ASSERT(n1bloc .eq. 1)
        ASSERT(n2bloc .eq. nueq)
!
        do j = 1, nueq
            hc = smdi(j)
            if (j .gt. 1) hc = hc-smdi(j-1)
            nbj = j-hc+1
            ASSERT(nbj .eq. 1)
!
! --------- MISE A ZERO PARTIE IMAGINAIRE DU MODE J
!
            do k = 1, neq
                vectass2(k) = dcmplx(vbasemo(1+(j-1)*neq+k-1), zero)
            end do
!
! --------- CALCUL PRODUIT MATRICE*MODE J
!
            call mcmult('ZERO', imatra, vectass2, vectass3, 1, &
                        .true._1)
            call zerlag(neq, zi(iddeeq), vectz=vectass3)
!
! --------- BOUCLE SUR LES INDICES DE LA COLONNE I
            do i = 1, nueq
!
! ------------ PRODUIT SCALAIRE VECTASS * MODE
                pij = dcmplx(zero, zero)
                do k = 1, neq
                    pij = pij+vectass3(k)*dcmplx(vbasemo(1+(i-1)*neq+k-1), zero)
                end do
!
! ------------ STOCKAGE DANS LE .UALF A LA BONNE PLACE (2 BLOCs)
                if (j .ge. i) then
                    zc(ldblo1+smdi(j)-j+i-1) = pij
                end if
                if (j .le. i) then
                    zc(ldblo2+smdi(i)-i+j-1) = pij
                end if
            end do
        end do
    end if
    AS_DEALLOCATE(vc=vectass2)
    AS_DEALLOCATE(vc=vectass3)
    AS_DEALLOCATE(vr=vbasemo)
!
!
    call ualfva(resu, 'G')
    call jedema()
end subroutine
