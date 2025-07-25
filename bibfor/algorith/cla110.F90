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

subroutine cla110(nomres, modgen)
    implicit none
!  P. RICHARD     DATE 22/04/91
!-----------------------------------------------------------------------
!  BUT : < MAILLAGE SQUELETTE SOUS-STRUCTURATION CLASSIQUE >
!  CREER LE MAILLAGE SQUELETTE CORRESPONDANT A UN MODELE GENERALISE
!
!-----------------------------------------------------------------------
!
! NOMRES  /I/ : NOM K8 DU MAILLAGE A CREER
! MODGEN  /I/ : NOM K8 DU MODELE GENERALISE
!
!
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/codlet.h"
#include "asterfort/compma.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvem.h"
#include "asterfort/getvtx.h"
#include "asterfort/gma110.h"
#include "asterfort/intet0.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/mgutdm.h"
#include "asterfort/nomcod.h"
#include "asterfort/pmppr.h"
#include "asterfort/r8inir.h"
#include "asterfort/recuma.h"
#include "asterfort/utmess.h"
#include "asterfort/uttrii.h"
#include "asterfort/wkvect.h"
!
!
!
!   PARAMETER : REPRESENTE LE NOMBRE MAX DE COMPOSANTES DE LA GRANDEUR
!   SOUS-JACENTE TRAITEE
!
    integer(kind=8) :: nbcmpm, ibid
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iad, iatyma, icomp, igd, igr, igrma
    integer(kind=8) :: ilstgr, ioc, iret, is, itcon, j, k
    integer(kind=8) :: l, ldcone, ldcoo, lddes, lddime, ldgrma
    integer(kind=8) :: ldskin, ldtitr, llcona, llma, llrot
    integer(kind=8) :: lltra, lltyp, lstac, ltdesc, ltfac, ltino, ltinv
    integer(kind=8) :: ltlima, ltlino, ltmail, ltnbgr, ltnbma, ltnbno, ltnogr
    integer(kind=8) :: ltnoma, ltrot, lttra, lutgma, lutnom, lutsst, maxgr
    integer(kind=8) :: maxma, nbcon, nbgr, nbgrut, nbincr, nbma, nbmat
    integer(kind=8) :: nbno, nbnot, nbskma, nbsst, nbstac, nbtemp, nbtgrm
    integer(kind=8) :: nbtmma, nbtmno, nbtmp, nbtout, nbuf, nbvgr, nbvma
    integer(kind=8) :: nctail, ngrma, ngrmat, ntail, nuact, numma, numno
    integer(kind=8) :: nusst
    real(kind=8) :: xnew
!-----------------------------------------------------------------------
    parameter(nbcmpm=10)
    character(len=8) :: nomres, modgen, tt, mailla, nomcou, nomsst
    character(len=16) :: css, cma, cgr, maicon, nomcon
    real(kind=8) :: xanc(3)
    character(len=24) :: repnom, modrot, modtra, gpptnm
    character(len=24) :: valk(2), nomgr
    real(kind=8) :: matrot(nbcmpm, nbcmpm)
    real(kind=8) :: matbuf(nbcmpm, nbcmpm), mattmp(nbcmpm, nbcmpm)
    character(len=8) :: k8bid, exclu
    integer(kind=8), pointer :: nldtyp(:) => null()
    real(kind=8), pointer :: vale(:) => null()
!
!-----------------------------------------------------------------------
    data tt/'&&CLA110'/
    data css, cma, cgr/'SOUS_STRUC', 'MAILLE', 'GROUP_MA'/
!-----------------------------------------------------------------------
!
    call jemarq()
    repnom = modgen//'      .MODG.SSNO'
    call jelira(repnom, 'NOMMAX', nbsst)
!
!-----PRISE EN COMPTE DE LA PRESENCE DES SST DANS LE SQUELETTE----------
!
    call wkvect(tt//'.ACTIF', 'V V I', nbsst, ltfac)
    call getfac(css, ioc)
!
    do i = 1, ioc
        call getvtx(css, 'NOM', iocc=i, scal=nomsst, nbret=ibid)
        call jenonu(jexnom(repnom, nomsst), nusst)
        if (nusst .eq. 0) then
            valk(1) = nomsst
            call utmess('A', 'ALGORITH12_49', sk=valk(1))
        else
            zi(ltfac+nusst-1) = 1
        end if
    end do
!
    nbstac = 0
    do i = 1, nbsst
        nbstac = nbstac+zi(ltfac+i-1)
    end do
!
    if (nbstac .eq. 0) then
        call utmess('F', 'ALGORITH12_50')
    end if
!
!-----DEFINITION DES REPERTOIRES DE TRAVAIL-----------------------------
    call wkvect(tt//'.DESC', 'V V I', nbsst, ltdesc)
    call wkvect(tt//'.NB.MA', 'V V I', nbstac, ltnbma)
    call wkvect(tt//'.NB.GR', 'V V I', nbstac, ltnbgr)
    do i = 1, nbstac
        zi(ltnbma-1+i) = 0
        zi(ltnbgr-1+i) = 0
    end do
    call jecreo(tt//'.NOM.SST', 'V N K24')
    call jeecra(tt//'.NOM.SST', 'NOMMAX', nbstac, ' ')
    call jecrec(tt//'.LISTE.MA', 'V V I', 'NU', 'DISPERSE', 'VARIABLE', &
                nbstac)
    call jecrec(tt//'.LISTE.NO', 'V V I', 'NU', 'DISPERSE', 'VARIABLE', &
                nbstac)
!
! --- RECHERCHE DES NOMS DE GROUPES DE MAILLES UTILISATEUR ---
    call getvtx(' ', 'EXCLUSIF', scal=exclu, nbret=ibid)
    call getfac('NOM_GROUP_MA', nbgrut)
    if (nbgrut .gt. 0) then
        call wkvect(tt//'.UT.NOM', 'V V K24', nbgrut, lutnom)
        call wkvect(tt//'.UT.SST', 'V V K8', nbgrut, lutsst)
        call wkvect(tt//'.UT.GMA', 'V V K24', nbgrut, lutgma)
        do i = 1, nbgrut
            call getvtx('NOM_GROUP_MA', 'NOM', iocc=i, scal=zk24(lutnom-1+i), nbret=ibid)
            call getvtx('NOM_GROUP_MA', 'SOUS_STRUC', iocc=i, scal=zk8(lutsst-1+i), nbret=ibid)
            call getvem(mailla, 'GROUP_MA', 'NOM_GROUP_MA', 'GROUP_MA', i, &
                        1, zk24(lutgma-1+i), ibid)
!           --- RECHERCHE SI LA SOUS-STRUCTURE EXISTE ---
            is = 0
90          continue
            is = is+1
            if (is .le. nbsst) then
                if (zi(ltfac-1+is) .ne. 0) then
                    call jenuno(jexnum(repnom, is), nomsst)
                    if (nomsst .ne. zk8(lutsst-1+i)) goto 90
                else
                    goto 90
                end if
            else
                valk(1) = zk8(lutsst-1+i)
                valk(2) = k8bid
                call utmess('F', 'ALGORITH12_51', nk=2, valk=valk)
            end if
        end do
    else
        lutsst = 1
        lutnom = 1
        lutgma = 1
    end if
!
!  ECRITURE DES NOMS DES SST ACTIVES
    icomp = 0
    do i = 1, nbsst
        if (zi(ltfac-1+i) .ne. 0) then
            call jenuno(jexnum(repnom, i), nomsst)
            call jecroc(jexnom(tt//'.NOM.SST', nomsst))
            icomp = icomp+1
            zi(ltdesc-1+i) = icomp
        else
            zi(ltdesc-1+i) = 0
        end if
    end do
!
!
!-----DETERMINATION DES DIMENSIONS MAX DES LISTES UTILISATEUR-----------
    maxma = 0
    maxgr = 0
!
    do i = 1, ioc
        call getvtx(css, cma, iocc=i, nbval=0, nbret=nbvma)
        maxma = max(maxma, -nbvma)
        call getvtx(css, cgr, iocc=i, nbval=0, nbret=nbvgr)
        maxgr = max(maxgr, -nbvgr)
    end do
!
!-----ALLOCATION VECTEUR DE TRAVAIL-------------------------------------
!
    ltnoma = 1
    ltnogr = 1
    ltmail = 1
    if (maxma .ne. 0) then
        call wkvect(tt//'.NOM.MA', 'V V K8', maxma, ltnoma)
    end if
    if (maxgr .ne. 0) then
        call wkvect(tt//'.NOM.GR', 'V V K24', maxgr, ltnogr)
    end if
    call wkvect(tt//'.MAILLAGE', 'V V K8', nbstac, ltmail)
!
!-----DETERMINATION DU NOMBRE DE MAILLES POUR CHAQUE SST ACTIVE---------
    ngrmat = 0
    do i = 1, ioc
        call getvtx(css, 'NOM', iocc=i, scal=nomsst, nbret=ibid)
        call jenonu(jexnom(repnom, nomsst), nusst)
        nuact = zi(ltdesc-1+nusst)
        call getvtx(css, cma, iocc=i, nbval=0, nbret=nbma)
        nbma = -nbma
        call getvtx(css, cgr, iocc=i, nbval=0, nbret=nbgr)
        nbgr = -nbgr
        call getvtx(css, 'TOUT', iocc=i, nbval=0, nbret=nbtout)
        nbtout = -nbtout
        call mgutdm(modgen, nomsst, ibid, 'NOM_MAILLAGE', ibid, &
                    mailla)
        zk8(ltmail+nuact-1) = mailla
        ngrma = 0
        if (nbtout .eq. 0) then
            call getvtx(css, cgr, iocc=i, nbval=nbgr, vect=zk24(ltnogr), &
                        nbret=ibid)
            call compma(mailla, nbgr, zk24(ltnogr), nbuf)
            nbma = nbma+nbuf
        else
            call dismoi('NB_MA_MAILLA', mailla, 'MAILLAGE', repi=nbma)
            call jeexin(mailla//'.GROUPEMA', iret)
            if (iret .ne. 0) then
                call jelira(mailla//'.GROUPEMA', 'NUTIOC', ngrma)
                call wkvect(tt//'.GR.'//nomsst, 'V V K24', ngrma, igrma)
            end if
            do igr = 1, ngrma
                call jenuno(jexnum(mailla//'.GROUPEMA', igr), nomgr)
                zk24(igrma-1+igr) = nomgr
            end do
        end if
        zi(ltnbma+nuact-1) = zi(ltnbma+nuact-1)+nbma
        zi(ltnbgr+nuact-1) = zi(ltnbgr+nuact-1)+ngrma
        ngrmat = ngrmat+ngrma
    end do
!
!
!-----ECRITURE ATTRIBUT LONGUEUR----------------------------------------
!
    do i = 1, nbstac
        ntail = zi(ltnbma+i-1)
        call jeecra(jexnum(tt//'.LISTE.MA', i), 'LONMAX', ntail, ' ')
        zi(ltnbma+i-1) = 0
    end do
!
!-----DETERMINATION DES LISTES DES MAILLES PAR SST ACTIVE---------------
!
    do i = 1, ioc
        call getvtx(css, 'NOM', iocc=i, scal=nomsst, nbret=ibid)
        call jenonu(jexnom(repnom, nomsst), nusst)
        nuact = zi(ltdesc-1+nusst)
        call getvtx(css, 'TOUT', iocc=i, nbval=0, nbret=nbtout)
        nbtout = -nbtout
        mailla = zk8(ltmail+nuact-1)
        call jeveuo(jexnum(tt//'.LISTE.MA', nuact), 'E', ltlima)
        if (nbtout .gt. 0) then
            call dismoi('NB_MA_MAILLA', mailla, 'MAILLAGE', repi=nbma)
            iad = ltlima+zi(ltnbma+nuact-1)
            do j = 1, nbma
                zi(iad+j-1) = j
            end do
            zi(ltnbma+nuact-1) = zi(ltnbma+nuact-1)+nbma
        else
            call getvtx(css, cma, iocc=i, nbval=0, nbret=nbma)
            nbma = -nbma
            call getvtx(css, cgr, iocc=i, nbval=0, nbret=nbgr)
            nbgr = -nbgr
            call getvtx(css, cma, iocc=i, nbval=nbma, vect=zk8(ltnoma), &
                        nbret=ibid)
            call getvtx(css, cgr, iocc=i, nbval=nbgr, vect=zk24(ltnogr), &
                        nbret=ibid)
            iad = ltlima+zi(ltnbma+nuact-1)
            call recuma(mailla, nbma, nbgr, zk8(ltnoma), zk24(ltnogr), &
                        nbskma, zi(iad))
            zi(ltnbma+nuact-1) = zi(ltnbma+nuact-1)+nbskma
        end if
!
        call jelibe(jexnum(tt//'.LISTE.MA', nuact))
    end do
    if (maxma .ne. 0) then
        call jedetr(tt//'.NOM.MA')
    end if
    if (maxgr .ne. 0) then
        call jedetr(tt//'.NOM.GR')
    end if
!
!-----TRI DES MAILLES ET COMPTAGE DES NOEUDS----------------------------
!
    call wkvect(tt//'.NB.NO', 'V V I', nbstac, ltnbno)
!
    do i = 1, nbstac
        call jeveuo(jexnum(tt//'.LISTE.MA', i), 'L', ltlima)
        nbtemp = zi(ltnbma+i-1)
        nbskma = nbtemp
        if (nbskma .ne. 0) call uttrii(zi(ltlima), nbskma)
        zi(ltnbma+i-1) = nbskma
        mailla = zk8(ltmail+i-1)
        maicon = mailla//'.CONNEX'
        icomp = 0
        do j = 1, nbskma
            numma = zi(ltlima+j-1)
            call jelira(jexnum(maicon, numma), 'LONMAX', nbno)
            icomp = icomp+nbno
        end do
        zi(ltnbno+i-1) = icomp
    end do
!
!-----ECRITURE ATTRIBUT DIMENSION DES NOEUDS----------------------------
!
    do i = 1, nbstac
        ntail = zi(ltnbno+i-1)
        call jeecra(jexnum(tt//'.LISTE.NO', i), 'LONMAX', ntail, ' ')
    end do
!
!-----RECUPERATION DES NOEUDS-------------------------------------------
!
    nbnot = 0
    nbmat = 0
    nctail = 0
!
    do i = 1, nbstac
        mailla = zk8(ltmail+i-1)
        maicon = mailla//'.CONNEX'
        call jeveuo(jexnum(tt//'.LISTE.MA', i), 'L', ltlima)
        call jeveuo(jexnum(tt//'.LISTE.NO', i), 'E', ltlino)
        nbma = zi(ltnbma+i-1)
        icomp = 0
        do j = 1, nbma
            numma = zi(ltlima-1+j)
            call jelira(jexnum(maicon, numma), 'LONMAX', nbtmp)
            call jeveuo(jexnum(maicon, numma), 'L', llma)
            nctail = nctail+nbtmp
            do k = 1, nbtmp
                icomp = icomp+1
                zi(ltlino+icomp-1) = zi(llma+k-1)
            end do
        end do
        call jelibe(maicon)
        call jelibe(jexnum(tt//'.LISTE.MA', i))
        nbtmp = icomp
        nbno = nbtmp
        if (nbno .ne. 0) call uttrii(zi(ltlino), nbno)
        call jelibe(jexnum(tt//'.LISTE.NO', i))
        zi(ltnbno+i-1) = nbno
        nbnot = nbnot+nbno
        nbmat = nbmat+nbma
    end do
!
!-----TRAITEMENT DES ORIENTATIONS ET DES TRANSLATIONS DES SST-----------
!
    call wkvect(tt//'.ROTATION', 'V V R', nbstac*3, ltrot)
    call wkvect(tt//'.TRANSLATION', 'V V R', nbstac*3, lttra)
    modrot = modgen//'      .MODG.SSOR'
    modtra = modgen//'      .MODG.SSTR'
    do i = 1, nbsst
        icomp = zi(ltdesc-1+i)
        if (icomp .ne. 0) then
            call jenuno(jexnum(repnom, i), nomsst)
            call jenonu(jexnom(modrot(1:19)//'.SSNO', nomsst), ibid)
            call jeveuo(jexnum(modrot, ibid), 'L', llrot)
            do k = 1, 3
                zr(ltrot+3*(icomp-1)+k-1) = zr(llrot+k-1)
            end do
        end if
    end do
    do i = 1, nbsst
        icomp = zi(ltdesc-1+i)
        if (icomp .ne. 0) then
            call jenuno(jexnum(repnom, i), nomsst)
            call jenonu(jexnom(modtra(1:19)//'.SSNO', nomsst), ibid)
            call jeveuo(jexnum(modtra, ibid), 'L', lltra)
            do k = 1, 3
                zr(lttra+3*(icomp-1)+k-1) = zr(lltra+k-1)
            end do
        end if
    end do
!
!-----ALLOCATION DES OBJETS MAILLAGE RESULTAT---------------------------
!
    nomcon = nomres//'.CONNEX'
    call wkvect(nomres//'.DIME', 'G V I', 6, lddime)
    call wkvect(nomres//'           .TITR', 'G V K80', 1, ldtitr)
    call wkvect(nomres//'         .NOMSST', 'G V K8', nbstac, lstac)
!
    call jecreo(nomres//'.NOMMAI', 'G N K8')
    call jeecra(nomres//'.NOMMAI', 'NOMMAX', nbmat)
    call jecrec(nomcon, 'G V I', 'NU', 'CONTIG', 'VARIABLE', &
                nbmat)
    call jeecra(nomcon, 'NUTIOC', nbmat)
    call jeecra(nomcon, 'LONT', nctail)
    call wkvect(nomres//'.TYPMAIL', 'G V I', nbmat, ibid)
!
    call jecreo(nomres//'.NOMNOE', 'G N K8')
    call jeecra(nomres//'.NOMNOE', 'NOMMAX', nbnot)
!
    gpptnm = nomres//'.PTRNOMMAI'
    call jecreo(gpptnm, 'G N K24')
    call jeecra(gpptnm, 'NOMMAX', nbstac+ngrmat)
    call jecrec(nomres//'.GROUPEMA', 'G V I', 'NO '//gpptnm, 'DISPERSE', 'VARIABLE', &
                nbstac+ngrmat)
!
    call wkvect(nomres//'.COORDO    .DESC', 'G V I', 3, lddes)
    call jeecra(nomres//'.COORDO    .DESC', 'DOCU', cval='CHGO')
    call wkvect(nomres//'.COORDO    .VALE', 'G V R', 3*nbnot, ldcoo)
!
!-----REMPLISSAGE DU TITRE----------------------------------------------
    zk80(ldtitr) = 'MAILLAGE SQUELETTE SOUS-STRUCTURATION CLASSIQUE'
!
!-----REMPLISSAGE DU DIME ET DU DESC------------------------------------
!
    call dismoi('NUM_GD', 'GEOM_R', 'GRANDEUR', repi=igd)
    zi(lddes) = igd
    zi(lddes+1) = -3
    zi(lddes+2) = 14
    zi(lddime) = nbnot
    zi(lddime+1) = 0
    zi(lddime+2) = nbmat
    zi(lddime+3) = 0
    zi(lddime+4) = 0
    zi(lddime+5) = 3
!
!-----ALLOCATION DU INV.SQUELETTE---------------------------------------
    call wkvect(nomres//'.INV.SKELETON', 'G V I', nbnot*2, ldskin)
!
!-----LET'S GET CRAZY !!!-----------------------------------------------
!
    nbtmma = 0
    nbtmno = 0
    nbtgrm = 0
    nbincr = 0
    itcon = 0
    call jeveuo(nomres//'.TYPMAIL', 'E', vi=nldtyp)
    call jeveuo(nomcon, 'E', ldcone)
!
    do i = 1, nbstac
        mailla = zk8(ltmail+i-1)
        maicon = mailla//'.CONNEX'
        call jeveuo(mailla//'.COORDO    .VALE', 'L', vr=vale)
        call jenuno(jexnum(tt//'.NOM.SST', i), nomsst)
        zk8(lstac-1+i) = nomsst
!
        call intet0(zr(ltrot+(i-1)*3), mattmp, 3)
        call intet0(zr(ltrot+(i-1)*3+1), matrot, 2)
        call r8inir(nbcmpm*nbcmpm, 0.d0, matbuf, 1)
        call pmppr(mattmp, nbcmpm, nbcmpm, 1, matrot, &
                   nbcmpm, nbcmpm, 1, matbuf, nbcmpm, &
                   nbcmpm)
        call r8inir(nbcmpm*nbcmpm, 0.d0, matrot, 1)
        call intet0(zr(ltrot+(i-1)*3+2), mattmp, 1)
        call pmppr(matbuf, nbcmpm, nbcmpm, 1, mattmp, &
                   nbcmpm, nbcmpm, 1, matrot, nbcmpm, &
                   nbcmpm)
!
        nbno = zi(ltnbno+i-1)
        call jeveuo(jexnum(tt//'.LISTE.NO', i), 'L', ltlino)
!
!  BOUCLE SUR LES NOEUDS GENERIQUES DE LA SST COURANTE
        call wkvect(tt//'.INV.MAILLA', 'V V I', nbnot, ltinv)
!
!-- L'ALLOCATION DOIT SE FAIRE POUR UNE LONGUEUR CORRESPONDANT
!-- AU PLUS GRAND NUMERO DE NOEUD
        ibid = 0
        do j = 1, nbno
            numno = zi(ltlino+j-1)
            if (numno .gt. ibid) ibid = numno
        end do
        call wkvect(tt//'.INV.NOEUD', 'V V I', ibid, ltino)
!
        do j = 1, nbno
            numno = zi(ltlino+j-1)
            nbtmno = nbtmno+1
            zi(ltinv-1+nbtmno) = numno
!-- NUMERO DES NOUVEAUX NOEUDS EN FONCTION DES ANCIENS
            zi(ltino+numno-1) = nbtmno
!
            if (nbtmno .le. 999999) then
                nomcou = 'N'
                call nomcod(nomcou, nbtmno, 2, 8)
            else
                nomcou = 'N'
                call codlet(nbtmno, 'D0', nomcou(2:8))
            end if
!
            call jecroc(jexnom(nomres//'.NOMNOE', nomcou))
            do k = 1, nbsst
                if (zi(ltdesc-1+k) .eq. i) zi(ldskin+nbtmno-1) = k
            end do
            zi(ldskin+nbnot+nbtmno-1) = numno
            do k = 1, 3
                xanc(k) = vale(1+(numno-1)*3+k-1)
            end do
            do k = 1, 3
                xnew = 0.d0
                do l = 1, 3
                    xnew = xnew+matrot(k, l)*xanc(l)
                end do
                zr(ldcoo+(nbtmno-1)*3+k-1) = xnew+zr(lttra+(i-1)*3+k-1)
            end do
        end do
!
        call jelibe(mailla//'.COORDO    .VALE')
        call jeveuo(jexnum(tt//'.LISTE.NO', i), 'L', ltlino)
!
!  BOUCLE SUR LES ELEMENTS GENERIQUES DE LA SST COURANTE
        nbma = zi(ltnbma+i-1)
        call jeveuo(jexnum(tt//'.LISTE.MA', i), 'L', ltlima)
        call jecroc(jexnom(nomres//'.GROUPEMA', nomsst))
        call jeecra(jexnom(nomres//'.GROUPEMA', nomsst), 'LONMAX', max(1, nbma))
        call jeecra(jexnom(nomres//'.GROUPEMA', nomsst), 'LONUTI', nbma)
        call jeveuo(jexnom(nomres//'.GROUPEMA', nomsst), 'E', ldgrma)
        nbtgrm = nbtgrm+1
        do j = 1, nbma
            numma = zi(ltlima+j-1)
            nbtmma = nbtmma+1
!
            if (nbtmma .le. 999999) then
                nomcou = 'MA'
                call nomcod(nomcou, nbtmma, 3, 8)
            else
                nomcou = 'M'
                call codlet(nbtmma, 'D0', nomcou(2:8))
            end if
!
            zi(ldgrma+j-1) = nbtmma
            call jecroc(jexnom(nomres//'.NOMMAI', nomcou))
            call jelira(jexnum(maicon, numma), 'LONMAX', nbcon)
            call jeecra(jexnum(nomcon, nbtmma), 'LONMAX', nbcon)
            call jeveuo(jexnum(maicon, numma), 'L', llcona)
!
            do k = 1, nbcon
                itcon = itcon+1
                zi(ldcone+itcon-1) = zi(ltino+zi(llcona+k-1)-1)
            end do
!
            call jeveuo(mailla//'.TYPMAIL', 'L', iatyma)
            lltyp = iatyma-1+numma
            nldtyp(nbtmma) = zi(lltyp)
        end do
!
        nbgr = zi(ltnbgr-1+i)
        if (nbgr .gt. 0) then
!       --- TRAITEMENT DES GROUPES DE MAILLES
            call jeveuo(tt//'.GR.'//nomsst, 'L', ilstgr)
            call gma110(nbgr, exclu, nbgrut, mailla, nomsst, &
                        nbtgrm, nomres, nbincr, zk24(ilstgr), zk8(lutsst), &
                        zk24(lutgma), zk24(lutnom))
        end if
        call jelibe(jexnum(tt//'.LISTE.MA', i))
        call jelibe(jexnom(nomres//'.GROUPEMA', nomsst))
        call jelibe(maicon)
        call jelibe(mailla//'.TYPMAIL')
        call jedetr(tt//'.INV.MAILLA')
        call jedetr(tt//'.INV.NOEUD')
!
        nbincr = nbincr+nbma
!
    end do
!
!
! --- MENAGE
!
    call jedetr(tt//'.ACTIF')
    call jedetr(tt//'.DESC')
    call jedetr(tt//'.NB.MA')
    call jedetr(tt//'.NB.GR')
    call jedetr(tt//'.NB.NO')
    call jedetr(tt//'.NOM.SST')
    call jedetr(tt//'.LISTE.MA')
    call jedetr(tt//'.LISTE.NO')
    call jedetr(tt//'.ROTATION')
    call jedetr(tt//'.TRANSLATION')
    call jedetr(tt//'.UT.NOM')
    call jedetr(tt//'.UT.SST')
    call jedetr(tt//'.UT.GMA')
    call jedetr(tt//'.GR.'//nomsst)
    call jedetr(tt//'.NOM.SST')
    call jedetr(tt//'.MAILLAGE')
!
    call jedema()
end subroutine
