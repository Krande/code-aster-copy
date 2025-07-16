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
subroutine irgmma(nomain, nomaou, nbmat, nummai, basz, &
                  nobj, nbel, versio)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/infniv.h"
#include "asterfort/irgmtb.h"
#include "asterfort/jecrec.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: nomain, nomaou
    integer(kind=8) :: nbmat, nummai(*), versio
    character(len=*) :: basz
!     TRANSFORME LE MAILLAGE "NOMAIN" EN UN MAILLAGE "NOMAOU"
!     LE MAILLAGE "NOMAOU" NE POSSEDE QUE DES MAILLES DE TYPE
!     POI1, SEG2, TRIA3, TETRA4 EN VERSION 1.0
!     + QUAD4, PENTA6, PYRAM5, HEXA8 EN VERSIO 1.2 (VOIR IRGMTB)
!     ------------------------------------------------------------------
    integer(kind=8) :: i, ima, nbma, nbmail, ifm, niv, ino, ima2, imav, iatyma
    integer(kind=8) :: jtitr
    integer(kind=8) :: jtypm, jdime, jopt, jnpt, nbmac, jmail, im, jnumol
    aster_logical :: logic
    character(len=1) :: base
    character(len=8) :: k8b, typm, typm2
    character(len=24) :: typmai, connex, nodime, cooval, coodsc
    character(len=24) :: titre, numold
    character(len=24) :: typmav, connev, nodimv, coovav, coodsv
    character(len=24) :: valk(2)
    integer(kind=8) :: ind, numel, nbcr, nbp
    integer(kind=8) :: nbmmax
    parameter(nbmmax=9999999)
!     --- TABLEAU DE DECOUPAGE
    integer(kind=8) :: ntyele, maxel, maxno
    parameter(ntyele=28)
    parameter(maxel=48)
    parameter(maxno=8)
    integer(kind=8), allocatable :: tdec(:, :, :)
    integer(kind=8) :: typd(ntyele, 3)
    integer(kind=8) :: vali(3)
!     NBRE, POINTEURS, NOM D'OBJET POUR CHAQUE TYPE D'ELEMENT
    integer(kind=8) :: nbel(ntyele), jel(ntyele), impr
    character(len=24) :: nobj(ntyele)
!     ------------------------------------------------------------------
!
    call infniv(ifm, niv)
    call jemarq()
    allocate (tdec(ntyele, maxel, maxno))
!
! --- INIT
    do i = 1, ntyele
        nbel(i) = 0
        jel(i) = 0
    end do
!
! --- TABLEAU DES INFOS DE DECOUPAGE
    call irgmtb(tdec, typd, versio)
!
    base = basz
!
    typmav = nomain//'.TYPMAIL        '
    connev = nomain//'.CONNEX         '
    nodimv = nomain//'.DIME           '
    coovav = nomain//'.COORDO    .VALE'
    coodsv = nomain//'.COORDO    .DESC'
!
    typmai = nomaou//'.TYPMAIL        '
    connex = nomaou//'.CONNEX         '
    nodime = nomaou//'.DIME           '
    cooval = nomaou//'.COORDO    .VALE'
    coodsc = nomaou//'.COORDO    .DESC'
    titre = nomaou//'           .TITR'
    numold = nomaou//'.NUMOLD         '
!
    call wkvect(titre, base//' V K80', 1, jtitr)
    zk80(jtitr) = 'MAILLAGE CREE PAR IRGMMA POUR GMSH'
!
    call jeveuo(typmav, 'L', jtypm)
    call jeveuo(nodimv, 'L', jdime)
    nbma = zi(jdime+3-1)
!
    logic = .false.
!
    if (nbmat .ne. 0) then
        nbmac = nbmat
        call wkvect('&&IRGMMA.NUME_MAILLE', 'V V I', nbmac, jmail)
        do ima = 1, nbmac
            zi(jmail+ima-1) = nummai(ima)
        end do
    else
        nbmac = nbma
        call wkvect('&&IRGMMA.NUME_MAILLE', 'V V I', nbmac, jmail)
        do ima = 1, nbmac
            zi(jmail+ima-1) = ima
        end do
    end if
!
! --- COMBIEN D'ELEMENTS DE CHAQUE TYPE VA-T-ON CREER ?
    do im = 1, nbmac
        ima = zi(jmail+im-1)
!
        ind = zi(jtypm+ima-1)
        call jenuno(jexnum('&CATA.TM.NOMTM', ind), typm)
!
! ---    NUMEL = EN QUOI ON DECOUPE, NBCR = COMBIEN ON EN CREER
        numel = typd(ind, 1)
        nbcr = typd(ind, 2)
        if (numel .ne. 0) then
            nbel(numel) = nbel(numel)+nbcr
        else
            call utmess('A', 'ALGELINE_64', sk=typm)
        end if
    end do
!
    nbmail = 0
    impr = 0
    do i = 1, ntyele
        nbmail = nbmail+nbel(i)
!
        if (nobj(i) .ne. ' ') then
            call wkvect(nobj(i), 'V V I', max(1, nbel(i)), jel(i))
            if (niv .ge. 1) then
                call jenuno(jexnum('&CATA.TM.NOMTM', i), typm)
                call jenuno(jexnum('&CATA.TM.NOMTM', typd(i, 1)), typm2)
                nbcr = typd(i, 2)
                nbp = typd(i, 3)
                if (nbel(i) .gt. 0) then
                    if (impr .eq. 0) then
                        call utmess('I', 'ALGELINE5_54')
                        impr = 1
                    end if
                    valk(1) = typm
                    valk(2) = typm2
                    vali(1) = nbel(i)
                    vali(2) = nbcr
                    vali(3) = nbp
                    call utmess('I', 'ALGELINE5_55', nk=2, valk=valk, ni=3, &
                                vali=vali)
                end if
            end if
        end if
    end do
!
    call wkvect(numold, 'V V I', max(1, nbmail), jnumol)
!
    call jedupo(nodimv, base, nodime, logic)
    call jedupo(coovav, base, cooval, logic)
    call jedupo(coodsv, base, coodsc, logic)
!
    call jeveuo(nodime, 'E', jdime)
    zi(jdime+3-1) = nbmail
!
! ----------------------------------------------------------------------
!     LE '.CONNEX'
! ----------------------------------------------------------------------
!
    call wkvect(typmai, base//' V I', nbmail, iatyma)
!
    call jecrec(connex, base//' V I', 'NU', 'CONTIG', 'VARIABLE', &
                nbmail)
!#MC  1*NBMAIL NE SUFFIT PAS ?
    call jeecra(connex, 'LONT', ntyele*nbmail, ' ')
!
    do i = 1, ntyele
        nbel(i) = 0
    end do
    imav = 0
!
    do im = 1, nbmac
        ima = zi(jmail+im-1)
!
        ind = zi(jtypm+ima-1)
        call jenuno(jexnum('&CATA.TM.NOMTM', ind), typm)
        call jeveuo(jexnum(connev, ima), 'L', jopt)
!
! ---    NUMEL = EN QUOI ON DECOUPE, NBCR = COMBIEN ON EN CREER
!        NBP = NBRE DE POINTS PAR ELEMENTS CREES
        numel = typd(ind, 1)
        nbcr = typd(ind, 2)
        nbp = typd(ind, 3)
        do i = 1, nbcr
            imav = imav+1
            if (imav .gt. nbmmax) then
                call codent(nbmmax, 'G', k8b)
                call utmess('F', 'ALGELINE_65', sk=k8b)
            end if
!
            ima2 = imav
            zi(iatyma-1+ima2) = numel
!    STOCKAGE DU NUMERO DE LA MAILLE INITIALE DANS NUMOLD POUR IRGMCE
            zi(jnumol-1+ima2) = ima
!
            call jeecra(jexnum(connex, ima2), 'LONMAX', nbp)
            call jeveuo(jexnum(connex, ima2), 'E', jnpt)
            do ino = 1, nbp
                zi(jnpt-1+ino) = zi(jopt-1+tdec(ind, i, ino))
            end do
            nbel(numel) = nbel(numel)+1
            zi(jel(numel)-1+nbel(numel)) = imav
        end do
!
    end do
!
    call jedetr('&&IRGMMA.NUME_MAILLE')
    deallocate (tdec)
!
    call jedema()
!
end subroutine
