! --------------------------------------------------------------------
! Copyright (C) 1991 - 201662025 - EDF R&D - www.code-aster.org
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

subroutine calprc(nomres, classe, basmod, nommat)
    implicit none
#include "jeveux.h"
#include "asterfort/copmod.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jexnum.h"
#include "asterfort/mcmult.h"
#include "asterfort/mtdscr.h"
#include "asterfort/mtexis.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/zeclag.h"
    character(len=24) :: nomres
    character(len=1) :: classe
    character(len=8) :: basmod
    character(len=19) :: nommat
!  BUT : < PROJECTION MATRICE SUR BASE QUELCONQUE >
!
!        CONSISTE A PROJETER UNE MATRICE ASSSEMBLEE COMPLEXE
!        SUR UNE BASE QUELCONQUE (PAS DE PROPRIETE D'ORTHOGONALITE)
!
!        LA MATRICE RESULTAT EST SYMETRIQUE ET STOCKEE TRIANGLE SUP
!
!-----------------------------------------------------------------------
!
! NOMRES /O/ : NOM K19 DE LA MATRICE CARREE RESULTAT
! CLASSE /I/ : CLASSE DE LA BASE JEVEUX DE L'OBJET RESULTAT
! BASMOD /I/ : NOM UTILISATEUR DE LA BASE MODALE DE PROJECTION
! NOMMAT /I/ : NOM UTITISATEUR DE LA MATRICE A PROJETER (RAIDEUR,MASSE)
!
!
!
!
    character(len=6) :: pgc
    character(len=14) :: num
    character(len=24) :: valk
    complex(kind=8) :: xprod, dcmplx, cbid
    integer(kind=8) :: ldref, nbdef, ntail, ldres, ier, lmat, neq, jdeb, ldres2
    integer(kind=8) :: iddeeq, idbase, ltvec1, ltvec2, i, j, k, iad, lddes, jrefa
    aster_logical :: lsym
    cbid = dcmplx(0.d0, 0.d0)
!
!
!-----------------------------------------------------------------------
    pgc = 'CALPRC'
!-----------------------------------------------------------------------
!
! --- CREATION DU .REFE
!
    call jemarq()
    call wkvect(nomres(1:18)//'_REFE', 'G V K24', 2, ldref)
    zk24(ldref) = basmod
    zk24(ldref+1) = nommat(1:8)
!
! --- RECUPERATION DES DIMENSIONS DE LA BASE MODALE
!
    call dismoi('NB_MODES_TOT', basmod, 'RESULTAT', repi=nbdef)
!
! --- ALLOCATION DE LA MATRICE RESULTAT
!
    call jeveuo(nommat(1:19)//'.REFA', 'L', jrefa)
    ntail = nbdef*(nbdef+1)/2
    if (zk24(jrefa-1+9) .eq. 'MS') then
        lsym = .true.
        call jecrec(nomres(1:18)//'_VALE', classe//' V C', 'NU', 'DISPERSE', &
                    'CONSTANT', 1)
    else
        lsym = .false.
        call jecrec(nomres(1:18)//'_VALE', classe//' V C', 'NU', 'DISPERSE', &
                    'CONSTANT', 2)
    end if
    call jeecra(nomres(1:18)//'_VALE', 'LONMAX', ntail)
    call jecroc(jexnum(nomres(1:18)//'_VALE', 1))
    call jeveuo(jexnum(nomres(1:18)//'_VALE', 1), 'E', ldres)
    if (.not. lsym) then
        call jecroc(jexnum(nomres(1:18)//'_VALE', 2))
        call jeveuo(jexnum(nomres(1:18)//'_VALE', 2), 'E', ldres2)
    end if
!
! --- CONTROLE D'EXISTENCE DE LA MATRICE
!
    call mtexis(nommat(1:8), ier)
    if (ier .eq. 0) then
        valk = nommat(1:8)
        call utmess('E', 'ALGORITH12_39', sk=valk)
    end if
!
! --- ALLOCATION DESCRIPTEUR DE LA MATRICE
!
    call mtdscr(nommat(1:8))
    call jeveuo(nommat(1:19)//'.&INT', 'E', lmat)
!
! --- RECUPERATION NUMEROTATION ET NB EQUATIONS
!
    call dismoi('NB_EQUA', nommat(1:8), 'MATR_ASSE', repi=neq)
    call dismoi('NOM_NUME_DDL', nommat(1:8), 'MATR_ASSE', repk=num)
    call jeveuo(num//'.NUME.DEEQ', 'L', iddeeq)
!
    call wkvect('&&'//pgc//'.BASEMO', 'V V R', nbdef*neq, idbase)
    call copmod(basmod, bmodr=zr(idbase), numer=num)
!
!
! --- ALLOCATION VECTEUR DE TRAVAIL
!
    call wkvect('&&'//pgc//'.VECT1', 'V V C', neq, ltvec1)
    call wkvect('&&'//pgc//'.VECT2', 'V V C', neq, ltvec2)
!
! --- PROJECTION SUR DEFORMEES
!
    do i = 1, nbdef
!
! ----- CALCUL PRODUIT MATRICE DEFORMEE
!
        do j = 1, neq
            zc(ltvec1+j-1) = dcmplx(zr(idbase+(i-1)*neq+j-1), 0.d0)
        end do
        call mcmult('ZERO', lmat, zc(ltvec1), zc(ltvec2), 1, &
                    .true._1)
        call zeclag(zc(ltvec2), neq, zi(iddeeq))
!        do j = 1, neq
!        end do
!
! ----- PRODUIT AVEC LA DEFORMEE COURANTE
!
        xprod = dcmplx(0.d0, 0.d0)
        do j = 1, neq
            xprod = xprod+zc(ltvec2-1+j)*dcmplx(zr(idbase+(i-1)*neq-1+ &
                                                   j), 0.d0)
        end do
!
        iad = i*(i+1)/2
        zc(ldres+iad-1) = xprod
        if (.not. lsym) zc(ldres2+iad-1) = xprod
        if (lsym) then
            jdeb = i+1
        else
            jdeb = 1
        end if
!
! ----- PRODUIT AVEC DEFORMEES D'ORDRE SUPERIEURE
!
!        if (i .lt. nbdef) then
        do j = jdeb, nbdef
            xprod = dcmplx(0.d0, 0.d0)
            do k = 1, neq
                xprod = xprod+zc(ltvec2-1+k)*dcmplx(zr(idbase+(j-1) &
                                                       *neq-1+k), 0.d0)
            end do
            if (j .gt. i) then
                iad = i+(j-1)*j/2
                zc(ldres+iad-1) = xprod
            else
                iad = j+(i-1)*i/2
                zc(ldres2+iad-1) = xprod
            end if
        end do
!        endif
!
    end do
!
    call jedetr('&&'//pgc//'.VECT1')
    call jedetr('&&'//pgc//'.VECT2')
    call jedetr('&&'//pgc//'.BASEMO')
!
! --- CREATION DU .DESC
!
    call wkvect(nomres(1:18)//'_DESC', 'G V I', 3, lddes)
    zi(lddes) = 2
    zi(lddes+1) = nbdef
    zi(lddes+2) = 2
!
    call jedema()
end subroutine
