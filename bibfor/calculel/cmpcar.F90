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
subroutine cmpcar(carte)
    implicit none
!
!     COMPRESSION D'1 CARTE :
! ( LORSQUE LES CMPS D'1 GRANDEUR N'ONT PAS ETE DONNEES SIMULTANEMENT)
!
!-----------------------------------------------------------------------
!
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/exisdg.h"
#include "asterfort/jacopo.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/meiden.h"
#include "asterfort/nbec.h"
#include "asterfort/scalai.h"
!
    character(len=19) :: carte
! ----------------------------------------------------------------------
!     ENTREES:
!       CARTE : NOM D'1 CARTE A COMPRIMER
!     SORTIES:
!      ON A RESTAURE LES OBJETS INITIAUX DE LA CARTE : .VALE,.NOLI,.LIMA
! ----------------------------------------------------------------------
!
!     VARIABLES LOCALES:
!     ------------------
    character(len=8) :: scal, ctype
    character(len=1) :: bas1
!
!
!     -- RECUPERATION DES OBJETS JEVEUX DE LA CARTE:
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i1, i2, i2lima, i3, i3desc, i3lima
    integer(kind=8) :: i3noli, i3vale, i4, iad, iad1, iad2
    integer(kind=8) :: iadgp, ianoma, iavale
    integer(kind=8) :: iavalp, iavtra, ico, icompt, iedit
    integer(kind=8) :: igd, ii, irtnu, isigne, j, k, n
    integer(kind=8) :: n1, nb, nbedi3, nbedit, nbmato, nboc, ncmp
    integer(kind=8) :: nec, num1, num2
    character(len=24), pointer :: noli(:) => null()
    integer(kind=8), pointer :: vret(:) => null()
    integer(kind=8), pointer :: lim2(:) => null()
    integer(kind=8), pointer :: numt(:) => null()
    integer(kind=8), pointer :: lipr(:) => null()
    integer(kind=8), pointer :: desc(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    call jeveuo(carte//'.DESC', 'L', vi=desc)
    call jeveuo(carte//'.VALE', 'L', iavale)
    call jeveuo(carte//'.VALP', 'L', iad1)
    call jelira(carte//'.VALP', 'TYPELONG', cval=ctype)
    call jeveuo(carte//'.NOMA', 'L', ianoma)
    call jelira(carte//'.DESC', 'CLAS', cval=bas1)
!
!
    igd = desc(1)
    nec = nbec(igd)
!     -- SCAL = I,R,C,K8,...
    scal = scalai(igd)
!
!     -- NCMP : NOMBRE MAXIMAL DE CMP POUR LA GRANDEUR.
!     ----------------------------------------------------
    call jelira(jexnum('&CATA.GD.NOMCMP', igd), 'LONMAX', ncmp)
!
!     -- RECUPERATION DES OBJETS  .NOLI .NUMT .VALP ET .DGP:
!     ------------------------------------------------------
    call jeveuo(carte//'.NOLI', 'L', vk24=noli)
    call jeveuo(carte//'.NUMT', 'L', vi=numt)
    call jeveuo(carte//'.VALP', 'L', iavalp)
    call jeveuo(carte//'.DGP ', 'L', iadgp)
    call jelira(carte//'.DGP ', 'LONMAX', n1)
!     NOMBRE TOTAL DE MAILLES:
!     CELLES DU MAILLAGE TOUTES LES SUPPL. DES LIGREL ATACHES A LA CARTE
    nbmato = n1/nec
!
!     -- ALLOCATION DES OBJETS: .LIPR .VRET .VTRA ET .LIM2 :
!     ------------------------------------------------------
    call jecreo(carte//'.LIPR', 'V V I')
    call jeecra(carte//'.LIPR', 'LONMAX', nbmato)
    call jeveuo(carte//'.LIPR', 'E', vi=lipr)
!     --LIPR CONTIENT A CHAQUE ITERARION, LA LISTE  DES
!     --MAILLES AFFECTEES A LA MEME GRANDEUR.
!
    call jecreo(carte//'.VRET', 'V V I')
    call jeecra(carte//'.VRET', 'LONMAX', nbmato)
    call jeveuo(carte//'.VRET', 'E', vi=vret)
!     --VRET NOTE POUR CHAQUE MAILLE SI ELLE A ETE RETENUE COMME MODELE
!
    call jecreo(carte//'.VTRA', 'V V L')
    call jeecra(carte//'.VTRA', 'LONMAX', nbmato)
    call jeveuo(carte//'.VTRA', 'E', iavtra)
    do i = 1, nbmato
        zl(iavtra-1+i) = .false.
    end do
!     --VTRA NOTE POUR CHAQUE MAILLE SI ELLE A ETE TRAITEE.
!
    call jecrec(carte//'.LIM2', 'V V I', 'NU', 'CONTIG', 'VARIABLE', &
                nbmato)
!     -- ON ESPERE QU'IL Y AURA MOINS DE GROUPES QUE DE MAILLES !!!
    call jeecra(carte//'.LIM2', 'LONT', nbmato, ' ')
    call jeveuo(carte//'.LIM2', 'E', vi=lim2)
!     --LIM2 EST LA COLLECTION QUI REMPLACERA .LIMA.
!
!
!     -- TRAITEMENT:
!     --------------
!
    nbedit = desc(3)
    ii = 0
    do iedit = 1, nbedit
        if (numt(3*(iedit-1)+3) .eq. 0) goto 11
        num1 = numt((iedit-1)*3+1)
        num2 = numt((iedit-1)*3+2)
        irtnu = 0
        do i = num1, num2
!           -- SI LA MAILLE A DEJA ETE TRAITEE: ON SORT DE LA BOUCLE.
            if (zl(iavtra-1+i)) goto 12
            icompt = 1
            vret(i) = iedit
            irtnu = irtnu+1
            lipr(icompt) = i
            i1 = iavalp-1+(i-1)*ncmp
            i2 = iadgp-1+(i-1)*nec
            do j = i+1, num2
                if (zl(iavtra-1+i)) goto 13
                i3 = iavalp-1+(j-1)*ncmp
                i4 = iadgp-1+(j-1)*nec
!              -- TESTE SI LES 2 GRANDEURS SONT PARFAITEMENT IDENTIQUES:
                if (meiden(scal(1:4), ncmp, i1, i3, nec, i2, i4)) then
                    icompt = icompt+1
                    lipr(icompt) = j
                    zl(iavtra-1+j) = .true.
                end if
13              continue
            end do
!           -- RECOPIE DE LA LISTE DE MAILLES .LIPR DANS .LIM2 :
!           -- ATTENTION .LIM2 CONTIENT LES NUMEROS TOTAUX DES MAILLES!
            call jecroc(jexnum(carte//'.LIM2', irtnu))
            call jeecra(jexnum(carte//'.LIM2', irtnu), 'LONMAX', icompt)
            do k = 1, icompt
                lim2(ii+k) = lipr(k)
            end do
            ii = ii+icompt
            zl(iavtra-1+i) = .true.
12          continue
        end do
11      continue
    end do
!
!     -- ON RECOPIE CE QU'IL FAUT DANS LES OBJETS FINAUX:
!     ---------------------------------------------------
!
    call jelira(carte//'.LIM2', 'NUTIOC', nbedi3)
    call jecreo(carte//'.DES3', 'V V I')
    call jeecra(carte//'.DES3', 'LONMAX', 3+nbedi3*(2+nec))
    call jeveuo(carte//'.DES3', 'E', i3desc)
    zi(i3desc-1+1) = desc(1)
    zi(i3desc-1+2) = nbedi3
    zi(i3desc-1+3) = nbedi3
!
    call jecreo(carte//'.NOL3', 'V V K24')
    call jeecra(carte//'.NOL3', 'LONMAX', nbedi3)
    call jeveuo(carte//'.NOL3', 'E', i3noli)
!
    call jecreo(carte//'.VAL3', 'V V '//scal(1:4))
    call jeecra(carte//'.VAL3', 'LONMAX', nbedi3*ncmp)
    call jeveuo(carte//'.VAL3', 'E', i3vale)
!
    call jecrec(carte//'.LIM3', 'V V I', 'NU', 'CONTIG', 'VARIABLE', &
                nbedi3)
    call jeecra(carte//'.LIM3', 'LONT', nbmato, ' ')
!
    icompt = 0
    do i = 1, nbmato
        if (vret(i) .le. 0) goto 21
        iedit = vret(i)
        icompt = icompt+1
!
!        --DES3 ET NOL3:
        zk24(i3noli-1+icompt) = zk24(i3noli-1+iedit)
        if (noli(iedit) (1:8) .eq. '        ') then
            zi(i3desc-1+3+2*(icompt-1)+1) = 3
        else
            zi(i3desc-1+3+2*(icompt-1)+1) = -3
        end if
        zi(i3desc-1+3+2*(icompt-1)+2) = icompt
        iad = i3desc-1+3+2*nbedi3+nec*(icompt-1)
        do k = 1, nec
            zi(iad+k) = zi(iadgp-1+(i-1)*nec+k)
        end do
!
!        --VAL3:
        ico = 0
        do k = 1, ncmp
            if (exisdg(zi(iadgp-1+(i-1)*nec+1), k)) then
                ico = ico+1
                call jacopo(1, ctype, iad1+ncmp*(i-1)+k-1, i3vale+ncmp*(icompt-1)+ico-1)
            end if
        end do
!
!        --LIMA:
        if (noli(iedit) (1:8) .eq. '        ') then
            isigne = 1
        else
            isigne = -1
        end if
        call jelira(jexnum(carte//'.LIM2', icompt), 'LONMAX', nb)
        call jeecra(jexnum(carte//'.LIM3', icompt), 'LONMAX', nb)
        call jeveuo(jexnum(carte//'.LIM2', icompt), 'L', i2lima)
        call jeveuo(jexnum(carte//'.LIM3', icompt), 'E', i3lima)
        num1 = numt((iedit-1)*3+1)
        do k = 1, nb
            zi(i3lima-1+k) = isigne*(zi(i2lima-1+k)-num1+1)
        end do
21      continue
    end do
!
!     ON DETRUIT LA CARTE INITIALE EST ON RECOPIE DEFINITIVEMENT:
!     -----------------------------------------------------------
!
    call jedetr(carte//'.DESC')
    call jedetr(carte//'.VALE')
    call jedetr(carte//'.NOLI')
    call jedetr(carte//'.LIMA')
!
!     DESC :
!     ------
    call jeveuo(carte//'.DES3', 'L', iad1)
    call jelira(carte//'.DES3', 'LONMAX', n)
    call jelira(carte//'.DES3', 'TYPELONG', cval=ctype)
    call jecreo(carte//'.DESC', bas1//' V I')
    call jeecra(carte//'.DESC', 'LONMAX', n)
    call jeecra(carte//'.DESC', 'DOCU', cval='CART')
    call jeveuo(carte//'.DESC', 'E', iad2)
    call jacopo(n, ctype, iad1, iad2)
!     NOLI :
!     ------
    call jelira(carte//'.NOL3', 'LONMAX', n)
    call jelira(carte//'.NOL3', 'TYPELONG', cval=ctype)
    call jeveuo(carte//'.NOL3', 'L', iad1)
    call jecreo(carte//'.NOLI', bas1//' V K24')
    call jeecra(carte//'.NOLI', 'LONMAX', n)
    call jeveuo(carte//'.NOLI', 'E', iad2)
    call jacopo(n, ctype, iad1, iad2)
!     VALE :
!     ------
    call jelira(carte//'.VAL3', 'LONMAX', n)
    call jelira(carte//'.VAL3', 'TYPELONG', cval=ctype)
    call jeveuo(carte//'.VAL3', 'E', iad1)
    call jecreo(carte//'.VALE', bas1//' V '//scal(1:4))
    call jeecra(carte//'.VALE', 'LONMAX', n)
    call jeveuo(carte//'.VALE', 'E', iad2)
    call jacopo(n, ctype, iad1, iad2)
!     LIMA :
!     ------
    call jelira(carte//'.LIM3', 'NMAXOC', nboc)
    call jelira(carte//'.LIM3', 'LONT', n)
    call jelira(carte//'.LIM3', 'TYPELONG', cval=ctype)
    call jeveuo(carte//'.LIM3', 'L', iad1)
    call jecrec(carte//'.LIMA', bas1//' V I', 'NU', 'CONTIG', 'VARIABLE', &
                nboc)
    call jeecra(carte//'.LIMA', 'LONT', n, ' ')
    call jeveuo(carte//'.LIMA', 'E', iad2)
!
    do i = 1, nboc
        call jelira(jexnum(carte//'.LIM3', i), 'LONMAX', nb)
        call jecroc(jexnum(carte//'.LIMA', i))
        call jeecra(jexnum(carte//'.LIMA', i), 'LONMAX', nb)
    end do
    call jacopo(n, ctype, iad1, iad2)
!
!        DESCRIPTION DE TOUS LES OBJETS DE TRAVAIL:
!
    call jedetr(carte//'.DES3')
    call jedetr(carte//'.DGP ')
    call jedetr(carte//'.LIM2')
    call jedetr(carte//'.LIM3')
    call jedetr(carte//'.LIPR')
    call jedetr(carte//'.NOL3')
    call jedetr(carte//'.NUMT')
    call jedetr(carte//'.VALP')
    call jedetr(carte//'.VAL3')
    call jedetr(carte//'.VRET')
    call jedetr(carte//'.VTRA')
    call jedema()
end subroutine
