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

subroutine cmqutr(basz, nomain, nomaou, nbma, nummai, &
                  prefix, ndinit)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/ingrma.h"
#include "asterfort/irgmtb.h"
#include "asterfort/jeccta.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/lxlgut.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: nbma, nummai(*), ndinit
    character(len=8) :: nomain, nomaou, prefix
    character(len=*) :: basz
!     OPTION = 'QUAD_TRIA3'
!
!     ------------------------------------------------------------------
    integer(kind=8) :: i, ima, nbmat, nbmail, typtri, nbtri, iret, nbgrno, nbnomx, nbpt
    integer(kind=8) :: ino, ima2, imav, iatyma, jvg, jtypm, jdime, jopt, jnpt
    integer(kind=8) :: nbno, ier, jgg, im, j, lgpref, nbmag, nbgrm, ifm, niv, iq4
    integer(kind=8) :: iq8, iq9, igrma, nbgm, jlgrma, jgrma, nbma2, jdec, ig, ind
    integer(kind=8) :: nbmais
    aster_logical :: logic
    character(len=1) :: base
    character(len=24) :: valk
    character(len=8) :: typm, nima
    character(len=24) :: typmai, connex, nodime, grpnoe, cooval
    character(len=24) :: coodsc, grpmai, nomg
    character(len=24) :: typmav, connev, nodimv, grpnov, gpptnn, coovav
    character(len=24) :: coodsv, grpmav, gpptnm
    integer(kind=8) :: versio
    parameter(versio=1)
!  --- TABLEAU DE DECOUPAGE
    integer(kind=8) :: ntyele, maxel, maxno
    parameter(ntyele=28)
    parameter(maxel=48)
    parameter(maxno=8)
    integer(kind=8), allocatable :: tdec(:, :, :)
    integer(kind=8) :: typd(ntyele, 3)
!     ------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
!
    allocate (tdec(ntyele, maxel, maxno))
!
!====
! 1. TABLEAU DE DECOUPAGE
!====
!
    call irgmtb(tdec, typd, versio)
!
!====
! 2. INITIALISATIONS DES NOMS D'OBJETS
!====
!
    base = basz
!
    typmav = nomain//'.TYPMAIL        '
    connev = nomain//'.CONNEX         '
    grpnov = nomain//'.GROUPENO       '
    grpmav = nomain//'.GROUPEMA       '
    nodimv = nomain//'.DIME           '
    coovav = nomain//'.COORDO    .VALE'
    coodsv = nomain//'.COORDO    .DESC'
!
    typmai = nomaou//'.TYPMAIL        '
    connex = nomaou//'.CONNEX         '
    grpnoe = nomaou//'.GROUPENO       '
    grpmai = nomaou//'.GROUPEMA       '
    nodime = nomaou//'.DIME           '
    cooval = nomaou//'.COORDO    .VALE'
    coodsc = nomaou//'.COORDO    .DESC'
!
    call jeveuo(typmav, 'L', jtypm)
    call jeveuo(nodimv, 'L', jdime)
!
!====
! 3. DIMENSIONNEMENT DU MAILLAGE RESULTAT
!    NBRE DE TRIANGLES A CREER
!====
!
!  NBMAT  : NB DE MAILLES DU MAILLAGE INITIAL
!  NBMA   : NB DE MAILLES POTENTIELLEMENT A DECOUPER
!  NBMAIL : NB DE MAILLES EN SORTIE DONT NBTRI TRIA3 CREES
    nbmat = zi(jdime+3-1)
!  --- VECTEUR A_DECOUPER_EN(NUM_MAILLE) = 0 OU N TRIA3 A CREER
    call wkvect('&&CMQUTR.A_DECOUPER_EN  ', 'V V I', nbmat, jdec)
!
    logic = .false.
    nbtri = 0
!
    iq4 = 0
    iq8 = 0
    iq9 = 0
    nbmail = nbmat
    do im = 1, nbma
        ima = nummai(im)
!
        call jenuno(jexnum('&CATA.TM.NOMTM', zi(jtypm+ima-1)), typm)
!
        if (typm .eq. 'QUAD4') then
            nbmail = nbmail-1
            nbtri = nbtri+2
            iq4 = iq4+1
            zi(jdec-1+ima) = 2
!
        else if (typm .eq. 'QUAD8') then
            nbmail = nbmail-1
            nbtri = nbtri+6
            iq8 = iq8+1
            zi(jdec-1+ima) = 6
!
        else if (typm .eq. 'QUAD9') then
            nbmail = nbmail-1
            nbtri = nbtri+6
            iq9 = iq9+1
            zi(jdec-1+ima) = 6
        end if
    end do
!
    if (niv .ge. 1) then
        write (ifm, 500) 1
        if (iq4 .ne. 0) write (ifm, 502) iq4, 'QUAD4', 2*iq4, 'TRIA3'
        if (iq8 .ne. 0) write (ifm, 502) iq8, 'QUAD8', 6*iq8, 'TRIA3'
        if (iq9 .ne. 0) write (ifm, 502) iq9, 'QUAD9', 6*iq9, 'TRIA3'
    end if
!
    nbmais = nbmail
    nbmail = nbmail+nbtri
!
    call jedupo(nodimv, base, nodime, logic)
    call jedupo(coovav, base, cooval, logic)
    call jedupo(coodsv, base, coodsc, logic)
!
!
    call jeveuo(nodime, 'E', jdime)
    zi(jdime+3-1) = nbmail
!
!====
! 4. CREATION DE SD DU MAILLAGE RESULTAT
!====
!
! 4.1. ==> CREATION DU .CONNEX
!
    call jenonu(jexnom('&CATA.TM.NOMTM', 'TRIA3'), typtri)
!
    call wkvect(typmai, base//' V I', nbmail, iatyma)
!
!     NBNOMX = NBRE DE NOEUDS MAX. POUR UNE MAILLE :
    call dismoi('NB_NO_MAX', '&CATA', 'CATALOGUE', repi=nbnomx)
!
    call jecrec(connex, base//' V I', 'NU', 'CONTIG', 'VARIABLE', &
                nbmail)
    call jeecra(connex, 'LONT', ival=nbnomx*nbmail)
!
! 4.2. ==> LE .GROUPMA EST CREE ICI,
!          LES GROUPES EUX-MEMES SERONT REMPLIS A LA VOLEE
!
    call jeexin(grpmav, igrma)
    if (igrma .ne. 0) then
        call jelira(grpmav, 'NOMUTI', ival=nbgrm)
        gpptnm = nomaou//'.PTRNOMMAI'
        call jecreo(gpptnm, 'G N K24')
        call jeecra(gpptnm, 'NOMMAX', ival=nbgrm)
        call jecrec(grpmai, base//' V I', 'NO '//gpptnm, 'DISPERSE', 'VARIABLE', &
                    nbgrm)
!     --- BCLE SUR LES GROUP_MA DU MAILLAGE INITIAL
        do i = 1, nbgrm
            call jenuno(jexnum(grpmav, i), nomg)
            call jeveuo(jexnum(grpmav, i), 'L', jgrma)
            call jelira(jexnum(grpmav, i), 'LONUTI', ival=nbmag)
            nbma2 = nbmag
!        --- BCLE SUR LES MAILLES DU GROUP_MA
            do j = 1, nbmag
                im = zi(jgrma-1+j)
                if (zi(jdec-1+im) .ne. 0) then
                    nbma2 = nbma2-1+zi(jdec-1+im)
                end if
            end do
            call jecroc(jexnom(grpmai, nomg))
!        --- LE NOUVEAU GROUP_MA CONTIENDRA NBMA2 MAILLES
            call jeecra(jexnom(grpmai, nomg), 'LONMAX', ival=max(1, nbma2))
            call jeecra(jexnom(grpmai, nomg), 'LONUTI', ival=nbma2)
            call jeecra(jexnom(grpmai, nomg), 'LONUTI', ival=0)
            if (niv .gt. 1) then
                write (ifm, *) 'GROUP_MA '//nomg, ' (', i, ') PASSE DE ', &
                    nbmag, ' A ', nbma2, ' MAILLES'
            end if
        end do
!     --- VECTEUR POUR STOCKER TEMPORAIREMENT LA LISTE DES GROUP_MA
!         D'UNE MAILLE
        call wkvect('&&CMQUTR.LISTE_GROUP_MA ', 'V V I', nbmag, jlgrma)
    end if
!
!====
! 5. ON PARCOURT LES MAILLES DU MAILLAGE INITIAL
!====
!
    lgpref = lxlgut(prefix)
    imav = 1
!
    do ima = 1, nbmat
!
        ind = zi(jtypm+ima-1)
        call jenuno(jexnum('&CATA.TM.NOMTM', ind), typm)
        call jeveuo(jexnum(connev, ima), 'L', jopt)
        call jelira(jexnum(connev, ima), 'LONMAX', ival=nbpt)
!
! 5.0. ==> PREPARE LA MISE DES GROUPES DE MAILLES
!
        nima = int_to_char8(ima)
        if (igrma .ne. 0) then
!        --- GROUP_MA CONTENANT IMA
            call ingrma(nomain, nima, zi(jlgrma), nbgm, ier)
        end if
!
! 5.1. ==> ON REGARDE SI LA MAILLE IMA DOIT ETRE DECOUPEE...
!
        if (zi(jdec-1+ima) .eq. 0) then
!
! 5.2. ==> ON CONSERVE LA MAILLE IMA TELLE QUELLE*
!          CAR IMA N'EST PAS DANS NUMMAI()
!
! 5.2.1. ==> TYPE DE MAILLE ET CONNECTIVITE
!
            ! ima2 = char8_to_int(nima)
            ima2 = imav
            imav = imav+1
            zi(iatyma-1+ima2) = zi(jtypm+ima-1)
!
            call jeecra(jexnum(connex, ima2), 'LONMAX', ival=nbpt)
            call jeveuo(jexnum(connex, ima2), 'E', jnpt)
            do ino = 1, nbpt
                zi(jnpt-1+ino) = zi(jopt+ino-1)
            end do
!
! 5.2.2. ==> MISE DES GROUPES DE MAILLES
!
            if (igrma .ne. 0 .and. ier .eq. 0 .and. nbgm .gt. 0) then
                do i = 1, nbgm
                    ig = zi(jlgrma-1+i)
                    call jeveuo(jexnum(grpmai, ig), 'E', jgrma)
                    call jelira(jexnum(grpmai, ig), 'LONUTI', ival=im)
                    im = im+1
!                  print *,'GROUP_MA ',IG,' : ',IM,' MAILLES'
                    zi(jgrma-1+im) = ima2
                    call jeecra(jexnum(grpmai, ig), 'LONUTI', ival=im)
                end do
            end if
!
! 5.3. ==> LA MAILLE IMA DOIT ETRE DECOUPE
!
        else
!
            nbpt = 3
            nbtri = zi(jdec-1+ima)
            do i = 1, nbtri
                ima2 = imav
                imav = imav+1
                zi(iatyma-1+ima2) = typtri
!
                call jeecra(jexnum(connex, ima2), 'LONMAX', ival=nbpt)
                call jeveuo(jexnum(connex, ima2), 'E', jnpt)
                do ino = 1, nbpt
!              --- TABLEAU DE DECOUPAGE SELON LE TYPE
                    zi(jnpt-1+ino) = zi(jopt-1+tdec(ind, i, ino))
                end do
!
                if (igrma .ne. 0 .and. ier .eq. 0 .and. nbgm .gt. 0) then
                    do j = 1, nbgm
                        ig = zi(jlgrma-1+j)
                        call jeveuo(jexnum(grpmai, ig), 'E', jgrma)
                        call jelira(jexnum(grpmai, ig), 'LONUTI', ival=im)
                        im = im+1
!                     print *,'GROUP_MA ',IG,' : ',IM,' MAILLES'
                        zi(jgrma-1+im) = ima2
                        call jeecra(jexnum(grpmai, ig), 'LONUTI', ival=im)
                    end do
                end if
!
            end do
!
        end if
!
!  --- MAILLE SUIVANTE
!
    end do
!
!====
! 6. LE .GROUPENO REPRIS A L'IDENTIQUE
!====
!
    call jeexin(grpnov, iret)
    if (iret .ne. 0) then
        call jelira(grpnov, 'NOMUTI', ival=nbgrno)
        gpptnn = nomaou//'.PTRNOMNOE'
        call jedetr(gpptnn)
        call jecreo(gpptnn, base//' N K24')
        call jeecra(gpptnn, 'NOMMAX', ival=nbgrno)
        call jecrec(grpnoe, base//' V I', 'NO '//gpptnn, 'DISPERSE', 'VARIABLE', &
                    nbgrno)
        do i = 1, nbgrno
            call jenuno(jexnum(grpnov, i), nomg)
            call jeveuo(jexnum(grpnov, i), 'L', jvg)
            call jelira(jexnum(grpnov, i), 'LONUTI', ival=nbno)
            call jeexin(jexnom(grpnoe, nomg), iret)
            if (iret .eq. 0) then
                call jecroc(jexnom(grpnoe, nomg))
            else
!           --- NE DEVRAIT PAS ARRIVER !
                valk = nomg
                call utmess('F', 'ALGELINE4_11', sk=valk)
            end if
            call jeecra(jexnom(grpnoe, nomg), 'LONMAX', ival=max(1, nbno))
            call jeecra(jexnom(grpnoe, nomg), 'LONUTI', ival=nbno)
            call jeveuo(jexnom(grpnoe, nomg), 'E', jgg)
            do j = 1, nbno
                zi(jgg-1+j) = zi(jvg-1+j)
            end do
        end do
    end if
!
!
!     -- RETASSAGE  DE CONNEX (QUI A ETE ALLOUEE TROP GRANDE) :
    call jeccta(connex)
    deallocate (tdec)
!
!
!
500 format('MOT CLE FACTEUR "MODI_MAILLE", OCCURRENCE ', i4)
502 format('  MODIFICATION DE ', i6, ' MAILLES ', a8,&
 &                     ' EN ', i6, ' MAILLES ', a8)
!
    call jedema()
!
end subroutine
