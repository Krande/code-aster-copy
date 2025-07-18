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

subroutine cmcrea(main, maout, nbocc)
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/cmfiss.h"
#include "asterfort/codent.h"
#include "asterfort/copisd.h"
#include "asterfort/cpclma.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
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
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbocc
    character(len=8) :: main, maout
!
! ----------------------------------------------------------------------
!  CREATION DE NOUVEAU MAILLAGE
! ----------------------------------------------------------------------
! IN        MAIN   K8  NOM DU MAILLAGE INITIAL
! IN/JXOUT  MAOUT  K8  NOM DU MAILLAGE TRANSFORME
! IN        NBOCC   I  NOMBRE D'OCCURENCES DES MOTS-CLES FACTEURS
! ----------------------------------------------------------------------
!
!
    character(len=16), parameter :: keywfact = 'CREA_FISS'
    integer(kind=8) :: nbnomx, nbmain, nbgmin, nbmaaj, nbgmaj, nbmato, nbgmto
    integer(kind=8) ::  numaco, numa, nbma, nbgm, nbno, ma, no, gm
    integer(kind=8) :: jadin, jadout
    integer(kind=8) :: jlnoma, jlconn, jltyma, jlngma, jlgpma
    integer(kind=8) :: jnoma, jtyma, jconn, jngma, jgpma
    integer(kind=8) :: jdim, ityin, ityout
    integer(kind=8) :: iret, iOcc, ib
    aster_logical :: false
!
    character(len=8) :: knum8
    character(len=24) :: linoma, liconn, lityma, lingma, ligpma
    character(len=24) :: valk, gno1, gno2
    character(len=24) :: dimin, dimout, tmain, tmaout, conin
    character(len=24) :: conout
    character(len=24) :: gmain, gmaout, nomgma, gmaptr
!
    data linoma/'&&CMCREA.LINOMA'/
    data liconn/'&&CMCREA.LICONN'/
    data lityma/'&&CMCREA.LITYMA'/
    data lingma/'&&CMCREA.LINGMA'/
    data ligpma/'&&CMCREA.LIGPMA'/
!
! NOMA  NOMS DES MAILLES CREEES (VECTEUR DE K8)
! CONN  CONNECTIVITE DES MAILLES CREEES (LISTE)
!          NBR DE NOEUDS DE LA MAILLE (ICI, TOUJOURS 4),
!          NUMEROS DES NDS DE LA MAILLE
! TYMA  TYPE DES MAILLES CREEES (VECTEUR I)
! NGMA  NOMS DES GFROUP_MA CREES
! GPMA  LISTE DES MAILLES DES GROUP_MA CREES
!          NBR DE MAILLES DU GROUP_MA
!          NUMERO DES MAILLES (NEGATIF QUAND NOUVELLE MAILLE)
! ----------------------------------------------------------------------
!
    call jemarq()
    false = .false.
!
! ----------------------------------------------------------------------
!                          INITIALISATION
! ----------------------------------------------------------------------
!
!    NOMBRE DE NOEUDS MAX. POUR UNE MAILLE :
    call dismoi('NB_NO_MAX', '&CATA', 'CATALOGUE', repi=nbnomx)
!
!    NOMBRE DE MAILLES, DE GROUP_MA
    call jeveuo(main//'.DIME', 'L', jadin)
    nbmain = zi(jadin-1+3)
    call jelira(main//'.GROUPEMA', 'NOMUTI', nbgmin)
!
!
!    LISTE DES OBJETS CREES PAR CHAQUE OCCURENCE DES MOTS-CLES
    call wkvect(linoma, 'V V K24', nbocc, jlnoma)
    call wkvect(liconn, 'V V K24', nbocc, jlconn)
    call wkvect(lityma, 'V V K24', nbocc, jltyma)
    call wkvect(lingma, 'V V K24', nbocc, jlngma)
    call wkvect(ligpma, 'V V K24', nbocc, jlgpma)
    do iOcc = 1, nbocc
        call codent(iOcc, 'D0', knum8)
        zk24(jlnoma-1+iOcc) = linoma(1:15)//'.'//knum8
        zk24(jlconn-1+iOcc) = liconn(1:15)//'.'//knum8
        zk24(jltyma-1+iOcc) = lityma(1:15)//'.'//knum8
        zk24(jlngma-1+iOcc) = lingma(1:15)//'.'//knum8
        zk24(jlgpma-1+iOcc) = ligpma(1:15)//'.'//knum8
    end do
!
!
! ----------------------------------------------------------------------
!    CREATION DES MAILLES, DES GROUP_MA, DES NOEUDS, DES GROUP_NO)
!           PARCOURS DES OCCURENCES DES MOTS-CLES FACTEURS
! ----------------------------------------------------------------------
!
    do iOcc = 1, nbocc
        call getvtx(keywfact, 'GROUP_NO_1', iocc=iOcc, scal=gno1, nbret=ib)
        call getvtx(keywfact, 'GROUP_NO_2', iocc=iOcc, scal=gno2, nbret=ib)
        call getvtx(keywfact, 'NOM', iocc=iOcc, scal=nomgma, nbret=ib)
!
        call cmfiss(main, gno1, gno2, nomgma, zk24(jlnoma-1+iOcc), &
                    zk24(jlconn-1+iOcc), zk24(jltyma-1+iOcc), &
                    zk24(jlngma-1+iOcc), zk24(jlgpma-1+iOcc))
    end do
!
!
! ----------------------------------------------------------------------
!                    DIMENSIONS DU NOUVEAU MAILLAGE
! ----------------------------------------------------------------------
!
    nbmaaj = 0
    nbgmaj = 0
    do iOcc = 1, nbocc
!
!      NOMBRE DE MAILLES AJOUTEES
        call jeexin(zk24(jlnoma-1+iOcc), iret)
        if (iret .ne. 0) then
            call jelira(zk24(jlnoma-1+iOcc), 'LONMAX', nbma)
            nbmaaj = nbmaaj+nbma
        end if
!
!      NOMBRE DE GROUP_MA AJOUTES
        call jeexin(zk24(jlngma-1+iOcc), iret)
        if (iret .ne. 0) then
            call jelira(zk24(jlngma-1+iOcc), 'LONMAX', nbgm)
            nbgmaj = nbgmaj+nbgm
        end if
!
    end do
!
    nbmato = nbmain+nbmaaj
    nbgmto = nbgmin+nbgmaj
!
!
!
! ----------------------------------------------------------------------
!                   CREATION DU NOUVEAU MAILLAGE
! ----------------------------------------------------------------------
!
! - CREATION DES OBJETS ET DUPLICATION DE LA PARTIE COMMUNE
!
!    OBJET .DIME
    dimin = main//'.DIME'
    dimout = maout//'.DIME'
    call jedupo(dimin, 'G', dimout, false)
    call jeveuo(dimout, 'E', jdim)
    zi(jdim-1+3) = nbmato
!
!    OBJET .TYPMAIL
    tmain = main//'.TYPMAIL'
    tmaout = maout//'.TYPMAIL'
    call wkvect(tmaout, 'G V I', nbmato, ityout)
    call jeveuo(tmain, 'L', ityin)
    do ma = 1, nbmain
        zi(ityout-1+ma) = zi(ityin-1+ma)
    end do
!
!    OBJET .CONNEX
    conin = main//'.CONNEX'
    conout = maout//'.CONNEX'
    call jecrec(conout, 'G V I', 'NU', 'CONTIG', 'VARIABLE', &
                nbmato)
    call jeecra(conout, 'LONT', nbnomx*nbmato, ' ')
    do ma = 1, nbmain
        call jelira(jexnum(conin, ma), 'LONMAX', nbno)
        call jeecra(jexnum(conout, ma), 'LONMAX', nbno)
        call jeveuo(jexnum(conin, ma), 'L', jadin)
        call jeveuo(jexnum(conout, ma), 'E', jadout)
        do no = 0, nbno-1
            zi(jadout+no) = zi(jadin+no)
        end do
    end do
!
!    OBJET .GROUPMA
    gmain = main//'.GROUPEMA'
    gmaout = maout//'.GROUPEMA'
    gmaptr = maout//'.PTRNOMMAI'
    call jecreo(gmaptr, 'G N K24')
    call jeecra(gmaptr, 'NOMMAX', nbgmto)
    call jecrec(gmaout, 'G V I', 'NO '//gmaptr, 'DISPERSE', 'VARIABLE', &
                nbgmto)
    do gm = 1, nbgmin
        call jenuno(jexnum(gmain, gm), nomgma)
        call jecroc(jexnom(gmaout, nomgma))
        call jelira(jexnum(gmain, gm), 'LONUTI', nbma)
        call jeecra(jexnom(gmaout, nomgma), 'LONMAX', nbma)
        call jeecra(jexnom(gmaout, nomgma), 'LONUTI', nbma)
        call jeveuo(jexnum(gmain, gm), 'L', jadin)
        call jeveuo(jexnom(gmaout, nomgma), 'E', jadout)
        do ma = 0, nbma-1
            zi(jadout+ma) = zi(jadin+ma)
        end do
    end do
!
!    DUPLICATION A L'IDENTIQUE .GROUPENO, .COORDO
!    (TANT QUE D'AUTRES MOTS CLES NE SONT PAS TRAITES)
    call cpclma(main, maout, 'GROUPENO', 'G')
    call copisd('CHAMP_GD', 'G', main//'.COORDO', maout//'.COORDO')
!
!    DUPLICATION A L'IDENTIQUE DES AUTRES OBJETS NON TRAITES
    call jedupo(main//'.NOMACR', 'G', maout//'.NOMACR', false)
    call jedupo(main//'.PARA_R', 'G', maout//'.PARA_R', false)
    call jedupo(main//'.SUPMAIL', 'G', maout//'.SUPMAIL', false)
    call jedupo(main//'.TYPL', 'G', maout//'.TYPL', false)
    call jedupo(main//'.ABSC_CURV', 'G', maout//'.ABSC_CURV', false)
!
!
! - AJOUT DES NOUVELLES MAILLES, DES NOUVEAUX GROUP_MA
!
    numaco = nbmain
    do iOcc = 1, nbocc
!
!  -    AJOUT DE NOUVELLES MAILLES
        nbma = 0
        call jeexin(zk24(jlnoma-1+iOcc), iret)
        if (iret .ne. 0) then
            call jelira(zk24(jlnoma-1+iOcc), 'LONMAX', nbma)
            call jeveuo(zk24(jlnoma-1+iOcc), 'L', jnoma)
            call jeveuo(zk24(jltyma-1+iOcc), 'L', jtyma)
            call jeveuo(zk24(jlconn-1+iOcc), 'L', jconn)
            do ma = 1, nbma
!
!          INSERTION DANS LE .TYPMAIL
                zi(ityout-1+numaco+ma) = zi(jtyma-1+ma)
!
!          INSERTION DANS LE .CONNEX
                nbno = zi(jconn)
                call jeecra(jexnum(conout, ma+numaco), 'LONMAX', nbno)
                call jeveuo(jexnum(conout, ma+numaco), 'E', jadout)
                do no = 1, nbno
                    zi(jadout-1+no) = zi(jconn+no)
                end do
                jconn = jconn+1+nbno
            end do
        end if
!
!
!   -   AJOUT DE NOUVEAUX GROUP_MA
        call jeexin(zk24(jlngma-1+iOcc), iret)
        if (iret .ne. 0) then
            call jelira(zk24(jlngma-1+iOcc), 'LONMAX', nbgm)
            call jeveuo(zk24(jlngma-1+iOcc), 'L', jngma)
            call jeveuo(zk24(jlgpma-1+iOcc), 'L', jgpma)
            do gm = 1, nbgm
!
!          INSERTION DANS LE .GROUPEMA
                nomgma = zk24(jngma-1+gm)
                call jeexin(jexnom(gmaout, nomgma), iret)
                if (iret .eq. 0) then
                    call jecroc(jexnom(gmaout, nomgma))
                else
                    valk = nomgma
                    call utmess('F', 'ALGELINE4_9', sk=valk)
                end if
!
                nbma = zi(jgpma)
                call jeecra(jexnom(gmaout, nomgma), 'LONMAX', nbma)
                call jeecra(jexnom(gmaout, nomgma), 'LONUTI', nbma)
                call jeveuo(jexnom(gmaout, nomgma), 'E', jadout)
!
                do ma = 1, nbma
                    numa = zi(jgpma+ma)
                    if (numa .le. 0) numa = -numa+numaco
                    zi(jadout-1+ma) = numa
                end do
                jgpma = jgpma+1+nbma
!
            end do
        end if
!
        numaco = numaco+nbma
    end do
!
!     -- RETASSAGE  DE CONOUT (QUI A ETE ALLOUEE TROP GRANDE) :
    call jeccta(conout)
!
! --- MENAGE
    call jedetr('&&CMCREA.LINOMA')
    call jedetr('&&CMCREA.LICONN')
    call jedetr('&&CMCREA.LITYMA')
    call jedetr('&&CMCREA.LINGMA')
    call jedetr('&&CMCREA.LIGPMA')
!
    call jedema()
end subroutine
