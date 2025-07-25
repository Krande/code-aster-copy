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

subroutine cmcovo(main, maout, nbma, lima, prefno, &
                  prefma, inima, epais, plan, trans)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8rddg.h"
#include "asterfort/codent.h"
#include "asterfort/codree.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
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
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: inima, nbma, lima(nbma)
    character(len=8) :: main, maout, prefno, prefma, plan, trans
    real(kind=8) :: epais
!         EXTRUSION DU GROUP_MA SURF EN GROUP_MA VOL
! ----------------------------------------------------------------------
! IN        MAIN   K8  NOM DU MAILLAGE INITIAL
! IN/JXOUT  MAOUT  K8  NOM DU MAILLAGE TRANSFORME
! IN        NBMA    I  NOMBRE DE MAILLES A TRAITER
! IN        LIMA    I  NUMERO ET TYPE DES MAILLES A TRAITER
! IN        PREFNO K8  PREFIXE DU NOM DES NOEUDS CREES (EX: N, NO, ...)
! IN        PREFMA K8  PREFIXE DU NOM DES MAILLES CREES (EX: N, NO, ..)
! IN        INIMA   I  NUMERO INITIAL DES MAILLES ET NOEUDS CREES
! IN        EPAIS   R  EPAISSEUR D'EXTRUSON
! IN        PLAN    K8 PLAN 'SUP' 'MOY' 'INF'
! IN        TRANS   K8 CAS PLAN ='MOY' ON TRANSLATE EN PEAU INF OU SUP
! ----------------------------------------------------------------------
    integer(kind=8) :: jdime, jcoor, nbnin, nbmin, nbnot, nbgrno, ifm, niv
    integer(kind=8) :: jnorn, ima, n1, n2, n3, nnoaj, ic, i, ij, iq4, it3
    integer(kind=8) :: jlisma, jnbnum, ino, jnorm
    integer(kind=8) :: jtypm, numa, nbno, lgno, inov, jnewm
    integer(kind=8) :: iret, jvale, kvale, lgnu, lgpref, nbgrmv
    integer(kind=8) :: typhex, typpen, iatyma, nbnomx, imav, lgnd, nbgrmn
    integer(kind=8) :: jopt, nbpt, jnpt, nbnuma, n4, jdimo, j, jvg
    integer(kind=8) :: nbmai, jgg, nbmat, jno, ima2
    character(len=8) :: knume, cdim, typm, ma1, ma2
    character(len=10) :: kangl
    character(len=24) :: normno, nonuma, grpmai, grpmav
    character(len=24) :: valk(4)
    character(len=24) :: typmav, connev, grpnov, nodimv
    character(len=24) :: coovav, coodsv, typmai
    character(len=24) :: connex, grpnoe, nodime, cooval, coodsc
    character(len=24) :: lisma, newma, grpnno, grpnma, nomg
    real(kind=8) :: coon1(3), coon2(3), coon3(3), coon4(3), n1n3(3), n1n2(3)
    real(kind=8) :: nx, ny, nz, nt(3), eps2, sinvec, cosvec
    real(kind=8) :: n4n2(3), n4n3(3), nq(3), norme, angl
    aster_logical :: logic
    character(len=24), pointer :: new_noeuds(:) => null()
    integer(kind=8), pointer :: noeuds(:) => null()
! ----------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
    logic = .false.
!
    typmav = main//'.TYPMAIL        '
    connev = main//'.CONNEX         '
    grpnov = main//'.GROUPENO       '
    grpmav = main//'.GROUPEMA       '
    nodimv = main//'.DIME           '
    coovav = main//'.COORDO    .VALE'
    coodsv = main//'.COORDO    .DESC'
!
    typmai = maout//'.TYPMAIL        '
    connex = maout//'.CONNEX         '
    grpnno = maout//'.PTRNOMNOE      '
    grpnma = maout//'.PTRNOMMAI      '
    grpnoe = maout//'.GROUPENO       '
    grpmai = maout//'.GROUPEMA       '
    nodime = maout//'.DIME           '
    cooval = maout//'.COORDO    .VALE'
    coodsc = maout//'.COORDO    .DESC'
!
    call jeveuo(typmav, 'L', jtypm)
!
    call jeveuo(nodimv, 'L', jdime)
    nbnin = zi(jdime)
    nbmin = zi(jdime+3-1)
!
!
    eps2 = epais/2.0d0
!
! --- RECUPERATION DU TABLEAU DES COORDONNEES :
!     ---------------------------------------
    call jeveuo(coovav, 'L', jcoor)
!
! --- RECUPERATION DU NOMBRE DE NOEUDS A AJOUTER
!     POUR DIMENSIONNER LES VECTEURS
!     -------------------------------------------
    AS_ALLOCATE(vk24=new_noeuds, size=nbnin)
    AS_ALLOCATE(vi=noeuds, size=nbnin)
    do ima = 1, nbma
        numa = lima(ima)
        call jenuno(jexnum('&CATA.TM.NOMTM', zi(jtypm+numa-1)), typm)
        if (typm .eq. 'QUAD4') then
        else if (typm .eq. 'TRIA3') then
        else
            call utmess('F', 'ALGELINE_14', sk=typm)
        end if
        call jelira(jexnum(connev, numa), 'LONMAX', nbno)
        call jeveuo(jexnum(connev, numa), 'L', jopt)
        do ino = 1, nbno
            noeuds(1+zi(jopt+ino-1)-1) = 1
        end do
    end do
!
    nnoaj = 0
    do ino = 1, nbnin
        if (noeuds(ino) .eq. 1) nnoaj = nnoaj+1
    end do
!
! --- STOCKAGE DES ELEMENTS DE TRAVAIL
!     --------------------------------
    normno = '&&CMCOVO.NORM'
    lisma = '&&CMCOVO.LISTE_MAILLES'
    nonuma = '&&CMCOVO.NB_NUME_MAILLE'
    newma = '&&CMCOVO.NEW_MAILLE'
!
    call wkvect(normno, 'V V R', 3*nbmin, jnorn)
    call wkvect(lisma, 'V V I', 27*nbnin, jlisma)
    call wkvect(nonuma, 'V V I', nbnin, jnbnum)
    call wkvect(newma, 'V V I', nbma, jnewm)
!
! --- STOCKAGE DES NOEUDS A TRAITER
!     ------------------------------
!
    do ima = 1, nbma
        numa = lima(ima)
        call jelira(jexnum(connev, numa), 'LONMAX', nbno)
        call jeveuo(jexnum(connev, numa), 'L', jopt)
!
! --- CALCUL DU PRODUIT VECTORIEL RELATIF A LA MAILLE NUMA
!     ----------------------------------------------------
        n1 = zi(jopt+1-1)
        n2 = zi(jopt+2-1)
        n3 = zi(jopt+3-1)
        if (nbno .eq. 4) n4 = zi(jopt+4-1)
        do ic = 1, 3
            coon1(ic) = zr(jcoor+3*(n1-1)+ic-1)
            coon2(ic) = zr(jcoor+3*(n2-1)+ic-1)
            coon3(ic) = zr(jcoor+3*(n3-1)+ic-1)
            if (nbno .eq. 4) then
                coon4(ic) = zr(jcoor+3*(n4-1)+ic-1)
            else
                coon4(ic) = 0.0d0
            end if
!
        end do
        n1n3 = coon3-coon1
        n1n2 = coon2-coon1
        call provec(n1n2, n1n3, nt)
!
        call normev(nt, norme)
!
        if (nbno .eq. 4) then
! --- ELEMENT QUAD: ON VERIFIE QUE LE NOEUD 4 N'EST PAS GAUCHE
! --- ON CALCULE LE PDV ET ON MOYENNE AVEC CELUI CALCULER PRECEDEMMENT
!
            n4n2 = coon2-coon4
            n4n3 = coon3-coon4
!
            call provec(n4n2, n4n3, nq)
!
            call normev(nq, norme)
! LA NORMALE EST INDICEE PAR LE NUMERO DE LA MAILLE NUMA.
            zr(jnorn+3*(numa-1)+1-1) = (nq(1)+nt(1))/2
            zr(jnorn+3*(numa-1)+2-1) = (nq(2)+nt(2))/2
            zr(jnorn+3*(numa-1)+3-1) = (nq(3)+nt(3))/2
        else
            zr(jnorn+3*(numa-1)+1-1) = nt(1)
            zr(jnorn+3*(numa-1)+2-1) = nt(2)
            zr(jnorn+3*(numa-1)+3-1) = nt(3)
        end if
    end do
!
! --- RECUPERATION DU NOMBRE DE MAILLES ET DE LA
!     LISTE DES MAILLES COMMUNE POUR UN NOEUD DONNE
!     -------------------------
    do ima = 1, nbma
        numa = lima(ima)
        call jelira(jexnum(connev, numa), 'LONMAX', nbno)
        call jeveuo(jexnum(connev, numa), 'L', jopt)
        do ij = 1, nbno
            ino = zi(jopt+ij-1)
            zi(jnbnum+ino-1) = zi(jnbnum+ino-1)+1
            zi(jlisma-1+27*(ino-1)+zi(jnbnum+ino-1)) = numa
        end do
    end do
!
! --- ON MOYENNE LES NORMALES
!
    call wkvect('&&CMCOVO.NORM_NO', 'V V R8', 3*nbnin, jnorm)
    do ino = 1, nbnin
        if (noeuds(ino) .eq. 0) cycle
        nbnuma = zi(jnbnum+ino-1)
        numa = zi(jlisma-1+27*(ino-1)+1)
        ma1 = int_to_char8(numa)
        zr(jnorm+3*(ino-1)) = zr(jnorn+3*(numa-1))
        zr(jnorm+3*(ino-1)+1) = zr(jnorn+3*(numa-1)+1)
        zr(jnorm+3*(ino-1)+2) = zr(jnorn+3*(numa-1)+2)
!
! ------ ON VERIFIE QUE L'ANGLE FORME PAR LES NORMALES < 90 DEGRES
!
        do ima = 2, nbnuma
            numa = zi(jlisma-1+27*(ino-1)+ima)
            cosvec = zr( &
                     jnorm+3*(ino-1))*zr(jnorn+3*(numa-1))+zr(jnorm+3*(ino-1)+1)*zr(jnorn+3*(&
                     &numa-1)+1)+zr(jnorm+3*(ino-1)+2)*zr(jnorn+3*(numa-1)+2 &
                     )
            call provec(zr(jnorm+3*(ino-1)), zr(jnorn+3*(numa-1)), nt)
            sinvec = nt(1)*nt(1)+nt(2)*nt(2)+nt(3)*nt(3)
            sinvec = sqrt(sinvec)
            angl = r8rddg()*atan2(sinvec, cosvec)
            if (abs(angl) .gt. 90.0d0) then
                nomg = int_to_char8(ino)
                ma2 = int_to_char8(numa)
                call codree(abs(angl), 'E', kangl)
                valk(1) = nomg
                valk(2) = ma1
                valk(3) = ma2
                valk(4) = kangl
                call utmess('A', 'ALGELINE_15', nk=4, valk=valk)
            end if
        end do
!
        do ima = 2, nbnuma
            numa = zi(jlisma-1+27*(ino-1)+ima)
            zr(jnorm+3*(ino-1)) = zr(jnorm+3*(ino-1))+zr(jnorn+3*(numa-1))
            zr(jnorm+3*(ino-1)+1) = zr(jnorm+3*(ino-1)+1)+zr(jnorn+3*(numa-1)+1)
            zr(jnorm+3*(ino-1)+2) = zr(jnorm+3*(ino-1)+2)+zr(jnorn+3*(numa-1)+2)
        end do
!
        zr(jnorm+3*(ino-1)) = zr(jnorm+3*(ino-1))/nbnuma
        zr(jnorm+3*(ino-1)+1) = zr(jnorm+3*(ino-1)+1)/nbnuma
        zr(jnorm+3*(ino-1)+2) = zr(jnorm+3*(ino-1)+2)/nbnuma
        call normev(zr(jnorm+3*(ino-1)), norme)
    end do
!
!
    nbmat = nbmin+nbma
    nbnot = nbnin+nnoaj
!
! ----------------------------------------------------------------------
!          ON AGRANDIT LE '.NOMNOE' ET LE '.COORDO    .VALE'
! ----------------------------------------------------------------------
!
    call jedupo(nodimv, 'G', nodime, logic)
    call jeveuo(nodime, 'E', jdimo)
    zi(jdimo) = nbnot
    zi(jdimo+2) = nbmat
    zi(jdimo+5) = 3
!
! --- TRAITEMENT DES NOEUDS AJOUTES
!
    lgno = lxlgut(prefno)
    inov = inima-1
    do ino = 1, nbnin
        if (noeuds(ino) .eq. 0) cycle
        inov = inov+1
        call codent(inov, 'G', knume)
        lgnu = lxlgut(knume)
!
        if (lgnu+lgno .gt. 8) then
            call utmess('F', 'ALGELINE_16')
        end if
    end do
!
! --- RECUPERATION DES COORDONNES DU MAIN ET CREATION
! --- DES COORDONNEES POUR MAOUT
!
    call jedupo(coodsv, 'G', coodsc, logic)
!
! --- MAOUT EST DE DIMENSION 3
!
    call jeveuo(coovav, 'L', jvale)
    call wkvect(cooval, 'G V R8', 3*nbnot, kvale)
    call codent(3, 'G', cdim)
    call jeecra(cooval, 'DOCU', cval=cdim)
!
! --- ON RECOPIE DANS MAOUT LES COORDONNEES DES NOEUDS DE MAIN
!
    do i = 1, 3*nbnin
        zr(kvale+i-1) = zr(jvale+i-1)
    end do
!
! --- POUR CHAQUE NOEUD A TRAITER ON TRANSLATE LES COORDONNEES
!     DU NOUVEAU NOEUDS DE N
!
    jno = nbnin
    do ino = 1, nbnin
        if (noeuds(ino) .eq. 0) cycle
        jno = jno+1
!
        nx = zr(jnorm+3*(ino-1))
        ny = zr(jnorm+3*(ino-1)+1)
        nz = zr(jnorm+3*(ino-1)+2)
!
        if (plan .eq. 'SUP') then
!
            zr(kvale+3*(jno-1)) = zr(kvale+3*(ino-1))-nx*epais
            zr(kvale+3*(jno-1)+1) = zr(kvale+3*(ino-1)+1)-ny*epais
            zr(kvale+3*(jno-1)+2) = zr(kvale+3*(ino-1)+2)-nz*epais
!
        else if (plan .eq. 'INF') then
!
            zr(kvale+3*(jno-1)) = zr(kvale+3*(ino-1))+nx*epais
            zr(kvale+3*(jno-1)+1) = zr(kvale+3*(ino-1)+1)+ny*epais
            zr(kvale+3*(jno-1)+2) = zr(kvale+3*(ino-1)+2)+nz*epais
!
        else if (plan .eq. 'MOY') then
!
            if (trans .eq. 'INF') then
                zr(kvale+3*(ino-1)) = zr(kvale+3*(ino-1))-nx*eps2
                zr(kvale+3*(ino-1)+1) = zr(kvale+3*(ino-1)+1)-ny*eps2
                zr(kvale+3*(ino-1)+2) = zr(kvale+3*(ino-1)+2)-nz*eps2
!
                zr(kvale+3*(jno-1)) = zr(kvale+3*(ino-1))+nx*epais
                zr(kvale+3*(jno-1)+1) = zr(kvale+3*(ino-1)+1)+ny*epais
                zr(kvale+3*(jno-1)+2) = zr(kvale+3*(ino-1)+2)+nz*epais
!
            else if (trans .eq. 'SUP') then
                zr(kvale+3*(ino-1)) = zr(kvale+3*(ino-1))+nx*eps2
                zr(kvale+3*(ino-1)+1) = zr(kvale+3*(ino-1)+1)+ny*eps2
                zr(kvale+3*(ino-1)+2) = zr(kvale+3*(ino-1)+2)+nz*eps2
!
                zr(kvale+3*(jno-1)) = zr(kvale+3*(ino-1))-nx*epais
                zr(kvale+3*(jno-1)+1) = zr(kvale+3*(ino-1)+1)-ny*epais
                zr(kvale+3*(jno-1)+2) = zr(kvale+3*(ino-1)+2)-nz*epais
!
            end if
        end if
    end do
!
! ----------------------------------------------------------------------
!     LE '.NOMMAI' ET LE '.CONNEX'
! ----------------------------------------------------------------------
    call jenonu(jexnom('&CATA.TM.NOMTM', 'HEXA8'), typhex)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'PENTA6'), typpen)
!
    call wkvect(typmai, 'G V I', nbmat, iatyma)
!
!     NBNOMX = NBRE DE NOEUDS MAX. POUR UNE MAILLE :
    call dismoi('NB_NO_MAX', '&CATA', 'CATALOGUE', repi=nbnomx)
!
    call jecrec(connex, 'G V I', 'NU', 'CONTIG', 'VARIABLE', &
                nbmat)
    call jeecra(connex, 'LONT', ival=nbnomx*nbmat)
!
! --- ON RECUPERE LES MAILLES DE MAIN DANS MAOUT
!
    do ima = 1, nbmin
!
        zi(iatyma-1+ima) = zi(jtypm+ima-1)
!
        call jeveuo(jexnum(connev, ima), 'L', jopt)
        call jelira(jexnum(connev, ima), 'LONMAX', nbpt)
!
        call jeecra(jexnum(connex, ima), 'LONMAX', nbpt)
        call jeveuo(jexnum(connex, ima), 'E', jnpt)
!
        do ino = 1, nbpt
            zi(jnpt-1+ino) = zi(jopt+ino-1)
        end do
    end do
!
! --- TRAITEMENT DES MAILLES AJOUTEES
!
    iq4 = 0
    it3 = 0
    lgpref = lxlgut(prefma)
    imav = inima-1
    do ima = 1, nbma
        numa = lima(ima)
        call jeveuo(jexnum(connev, numa), 'L', jopt)
        call jenuno(jexnum('&CATA.TM.NOMTM', zi(jtypm+numa-1)), typm)
        imav = imav+1
        call codent(imav, 'G', knume)
        lgnd = lxlgut(knume)
        if (lgnd+lgpref .gt. 8) then
            call utmess('F', 'ALGELINE_17')
        end if
        nomg = prefma(1:lgpref)//knume
!
        ima2 = char8_to_int(nomg)
        zi(jnewm+ima-1) = ima2
!
        if (typm .eq. 'QUAD4') then
!             -----------------
! --------- CREATION DE LA MAILLE HEXA8
            nbpt = 8
            zi(iatyma-1+ima2) = typhex
!
            call jeecra(jexnum(connex, ima2), 'LONMAX', nbpt)
            call jeveuo(jexnum(connex, ima2), 'E', jnpt)
            do ino = 1, 4
                zi(jnpt-1+ino) = zi(jopt-1+ino)
            end do
            do ino = 5, 8
                zi(jnpt-1+ino) = char8_to_int(new_noeuds(1+zi(jopt-1+ino-4)-1))
            end do
            iq4 = iq4+1
!
!
        else if (typm .eq. 'TRIA3') then
!             -----------------
! --------- CREATION DE LA MAILLE PENTA6
            nbpt = 6
            zi(iatyma-1+ima2) = typpen
!
            call jeecra(jexnum(connex, ima2), 'LONMAX', nbpt)
            call jeveuo(jexnum(connex, ima2), 'E', jnpt)
            do ino = 1, 3
                zi(jnpt-1+ino) = zi(jopt-1+ino)
            end do
            do ino = 4, 6
                zi(jnpt-1+ino) = char8_to_int(new_noeuds(1+zi(jopt-1+ino-3)-1))
!
            end do
            it3 = it3+1
!
        end if
!
    end do
!  -------------------------------------------------------------
!                       CREATION DES GROUP_MA
!  -------------------------------------------------------------
    call jeexin(grpmav, iret)
    if (iret .eq. 0) then
        nbgrmv = 0
    else
        call jelira(grpmav, 'NOMUTI', nbgrmv)
    end if
    nbgrmn = nbgrmv+1
    if (nbgrmn .ne. 0) then
        call jecreo(grpnma, 'G N K24')
        call jeecra(grpnma, 'NOMMAX', nbgrmn)
        call jecrec(grpmai, 'G V I', 'NO '//grpnma, 'DISPERSE', 'VARIABLE', &
                    nbgrmn)
        do i = 1, nbgrmv
            call jenuno(jexnum(grpmav, i), nomg)
            call jeexin(jexnom(grpmai, nomg), iret)
            if (iret .eq. 0) then
                call jecroc(jexnom(grpmai, nomg))
            else
                valk(1) = nomg
                call utmess('F', 'ALGELINE4_9', sk=valk(1))
            end if
            call jeveuo(jexnum(grpmav, i), 'L', jvg)
            call jelira(jexnum(grpmav, i), 'LONUTI', nbmai)
            call jeecra(jexnom(grpmai, nomg), 'LONMAX', max(1, nbmai))
            call jeecra(jexnom(grpmai, nomg), 'LONUTI', nbmai)
            call jeveuo(jexnom(grpmai, nomg), 'E', jgg)
            do j = 0, nbmai-1
                zi(jgg+j) = zi(jvg+j)
            end do
        end do
!
        call getvtx('COQU_VOLU', 'NOM', iocc=1, scal=nomg, nbret=n1)
        call jeexin(jexnom(grpmai, nomg), iret)
        if (iret .eq. 0) then
            call jecroc(jexnom(grpmai, nomg))
        else
            valk(1) = nomg
            call utmess('F', 'ALGELINE4_9', sk=valk(1))
        end if
        call jeecra(jexnom(grpmai, nomg), 'LONMAX', max(1, nbma))
        call jeecra(jexnom(grpmai, nomg), 'LONUTI', nbma)
        call jeveuo(jexnom(grpmai, nomg), 'E', jgg)
        do j = 1, nbma
            zi(jgg+j-1) = zi(jnewm+j-1)
        end do
    end if
!  -------------------------------------------------------------
!                       CREATION DES GROUP_NO
!  -------------------------------------------------------------
    call jeexin(grpnov, iret)
    if (iret .ne. 0) then
        call jelira(grpnov, 'NOMUTI', nbgrno)
        call jecreo(grpnno, 'G N K24')
        call jeecra(grpnno, 'NOMMAX', nbgrno)
        call jecrec(grpnoe, 'G V I', 'NO '//grpnno, 'DISPERSE', 'VARIABLE', &
                    nbgrno)
        do i = 1, nbgrno
            call jenuno(jexnum(grpnov, i), nomg)
            call jeveuo(jexnum(grpnov, i), 'L', jvg)
            call jelira(jexnum(grpnov, i), 'LONUTI', nbno)
            call jeexin(jexnom(grpnoe, nomg), iret)
            if (iret .eq. 0) then
                call jecroc(jexnom(grpnoe, nomg))
            else
                valk(1) = nomg
                call utmess('F', 'ALGELINE4_11', sk=valk(1))
            end if
            call jeecra(jexnom(grpnoe, nomg), 'LONMAX', max(1, nbno))
            call jeecra(jexnom(grpnoe, nomg), 'LONUTI', nbno)
            call jeveuo(jexnom(grpnoe, nomg), 'E', jgg)
            do j = 0, nbno-1
                zi(jgg+j) = zi(jvg+j)
            end do
        end do
    end if
!
!     -- RETASSAGE  DE CONNEX (QUI A ETE ALLOUEE TROP GRANDE) :
    call jeccta(connex)
!
! ----------------------------------------------------------------------
!
    if (niv .ge. 1) then
        write (ifm, 100) 1
        if (iq4 .ne. 0) write (ifm, 102) iq4
        if (it3 .ne. 0) write (ifm, 104) it3
    end if
!
    AS_DEALLOCATE(vk24=new_noeuds)
    AS_DEALLOCATE(vi=noeuds)
    call jedetr('&&CMCOVO.TRAV')
    call jedetr('&&CMCOVO.NORM_NO')
!     --------------------------------
    call jedetr(normno)
    call jedetr(lisma)
    call jedetr(nonuma)
    call jedetr(newma)
!
100 format('MOT CLE FACTEUR "COQU_VOLU", OCCURRENCE ', i4)
102 format('  EXTRUSION DE ', i6, ' MAILLES "QUAD4" EN "HEXA8"')
104 format('  EXTRUSION DE ', i6, ' MAILLES "TRIA3" EN "PENTA6"')
!
    call jedema()
!
end subroutine
