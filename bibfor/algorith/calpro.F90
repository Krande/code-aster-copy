! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine calpro(nomres, classe, basmod, nommat)
    implicit none
! P. RICHARD     DATE 23/05/91
!-----------------------------------------------------------------------
!
!  BUT : < PROJECTION MATRICE SUR BASE QUELCONQUE >
!
!        CONSISTE A PROJETER UNE MATRICE ASSSEMBLEE SUR UNE BASE
!        QUELCONQUE (PAS DE PROPRIETE D'ORTHOGONALITE)
!
!        LA MATRICE RESULTAT EST SYMETRIQUE ET STOCKEE TRIANGLE SUP
!
!-----------------------------------------------------------------------
!
! NOMRES /I/ : NOM K19 DE LA MATRICE CARREE RESULTAT
! CLASSE /I/ : CLASSE DE LA BASE JEVEUX DE L'OBJET RESULTAT
! BASMOD /I/ : NOM UT DE LA BASE MODALE DE PROJECTION
! NOMMAT /I/ : NOM UT DE LA MATRICE A PROJETER (RAIDEUR,MASSE)
!
!
!
!
#include "jeveux.h"
#include "asterfort/copmod.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jexnum.h"
#include "asterfort/mrmult.h"
#include "asterfort/mtdscr.h"
#include "asterfort/mtexis.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/zerlag.h"
#include "blas/ddot.h"
    character(len=1) :: classe, typ1
    character(len=6) :: pgc
    character(len=8) :: basmod
    character(len=19) :: nommat
    character(len=14) :: num
    character(len=24) :: nomres
    character(len=24) :: valk
    character(len=24) :: nomcha
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer :: i, iad, idbase,  ier, iret, jrefa, jdeb
    integer :: j, lddes, ldref, ldres, ldres2, lmat, ltvec1, nbdef
    integer :: neq, ntail
    aster_logical :: lsym
    real(kind=8) :: xprod
    complex(kind=8) :: cbid
    integer, pointer :: deeq(:) => null()
    cbid = dcmplx(0.d0, 0.d0)
!-----------------------------------------------------------------------
    pgc='CALPRO'
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
! ----VERIFICATION DU TYPE DES VECTEURS PROPRES DANS LA BASE
!
    call rsexch('F', basmod, 'DEPL', 1, nomcha,&
                iret)
    call jelira(nomcha(1:19)//'.VALE', 'TYPE', cval=typ1)
    if (typ1 .eq. 'C') then
        valk = basmod
        call utmess('F', 'ALGORITH12_16', sk=valk)
    endif
!
! --- ALLOCATION DE LA MATRICE RESULTAT
!
    call jeveuo(nommat(1:19)//'.REFA', 'L', jrefa)
    ntail = nbdef* (nbdef+1)/2
    if (zk24(jrefa-1+9) .eq. 'MS') then
       lsym = .true.
       call jecrec(nomres(1:18)//'_VALE', classe//' V R', 'NU', 'DISPERSE', &
                   'CONSTANT',1)
    else
       lsym = .false.
       call jecrec(nomres(1:18)//'_VALE', classe//' V R', 'NU', 'DISPERSE', &
                   'CONSTANT',2)
    endif
    call jeecra(nomres(1:18)//'_VALE', 'LONMAX', ntail)
    call jecroc(jexnum(nomres(1:18)//'_VALE', 1))
    call jeveuo(jexnum(nomres(1:18)//'_VALE', 1), 'E', ldres)
    if (.not.lsym) then
        call jecroc(jexnum(nomres(1:18)//'_VALE', 2))
        call jeveuo(jexnum(nomres(1:18)//'_VALE', 2), 'E', ldres2)
    endif
!
! --- CONTROLE D'EXISTENCE DE LA MATRICE
!
    call mtexis(nommat(1:8), ier)
    if (ier .eq. 0) then
        valk = nommat(1:8)
        call utmess('E', 'ALGORITH12_39', sk=valk)
    endif
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
    call jeveuo(num//'.NUME.DEEQ', 'L', vi=deeq)
!
    call wkvect('&&CALPRO.BASEMO', 'V V R', nbdef*neq, idbase)
    call copmod(basmod, numer=num, bmodr=zr(idbase))
!
!
! --- ALLOCATION VECTEUR DE TRAVAIL
!
    call wkvect('&&'//pgc//'.VECT1', 'V V R', neq, ltvec1)
!
! --- PROJECTION SUR DEFORMEES
!
    do i = 1, nbdef
!
! ----- CALCUL PRODUIT MATRICE DEFORMEE
!
        call mrmult('ZERO', lmat, zr(idbase+(i-1)*neq), zr(ltvec1), 1,&
                    .true._1)
        call zerlag(neq, deeq, vectr=zr(ltvec1))
!
! ----- PRODUIT AVEC LA DEFORMEE COURANTE
!
         xprod= ddot(neq,zr(ltvec1),1,zr(idbase+(i-1)*neq),1)
         iad = i*(i+1)/2
         zr(ldres+iad-1) = xprod
         if (.not.lsym) zr(ldres2+iad-1) = xprod
         if (lsym) then
          jdeb = i+1
         else
          jdeb=1
         endif
!
! ----- PRODUIT AVEC DEFORMEES D'ORDRE SUPERIEURE
!
!       if (i .lt. nbdef) then
            do j = jdeb, nbdef
                xprod= ddot(neq,zr(ltvec1),1,zr(idbase+(j-1)*neq),1)
                if (j.gt.i) then
                  iad = i+(j-1)*j/2
                  zr(ldres+iad-1) = xprod
                else
                  iad = j+(i-1)*i/2
                  zr(ldres2+iad-1) = xprod
                end if
            end do
!       endif
!
    end do
!
    call jedetr('&&'//pgc//'.VECT1')
!
! --- CREATION DU .DESC
!
    call wkvect(nomres(1:18)//'_DESC', 'G V I', 3, lddes)
    zi(lddes) = 2
    zi(lddes+1) = nbdef
    zi(lddes+2) = 2
!
!
! --- MENAGE
!
    call jedetr('&&CALPRO.BASEMO')
!
    call jedema()
end subroutine
