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
subroutine ircers(ifi, ligrel, nbgrel, longr, ncmpmx, &
                  vale, nomgd, nomcmp, titr, nomel, &
                  loc, celd, nbnoma, permut, maxnod, &
                  typma, nomsd, nomsym, ir, nbmat, &
                  nummai, lmasu, ncmpu, nucmp, nbcmp, &
                  ncmps, nocmpl)
! aslint: disable=W1504
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/dgmode.h"
#include "asterfort/digdel.h"
#include "asterfort/ecrtes.h"
#include "asterfort/exisdg.h"
#include "asterfort/irgags.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/lxlgut.h"
#include "asterfort/lxliis.h"
#include "asterfort/nbec.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: maxnod, ifi, ligrel(*), nbgrel, longr(*), ncmpmx, celd(*), ncmpu
    integer(kind=8) :: nucmp(*), nbnoma(*), typma(*), permut(maxnod, *), nbmat
    integer(kind=8) :: nummai(*), nbcmp, ncmps(*)
    character(len=*) :: nomgd, nomcmp(*), nomel(*), loc, titr, nomsym, nomsd
    character(len=*) :: nocmpl(*)
    real(kind=8) :: vale(*)
    aster_logical :: lmasu
!        ECRITURE D'UN CHAMELEM SUR FICHIER UNIVERSEL, DATASET TYPE 56
!                                                                OU 57
!        A VALEURS REELLES
!  ENTREE:
!     IFI   : UNITE LOGIQUE DU FICHIER UNIVERSEL
!     LIGREL: LIGREL COMPLET
!     NBGREL: NOMBRE DE GRELS
!     LONGR : POINTEUR DE LONGUEUR DE LIGREL
!     NCMPMX: NOMBRE MAXI DE CMP DE LA GRANDEUR NOMGD
!     VALE  : VALEURS DU CHAM_NO
!     NOMGD : NOM DE LA GRANDEUR: SIEF_R, EPSI_R,...
!     NOMCMP: NOMS DES CMP
!     TITR  : 1 LIGNE  DE TITRE
!     NOMEL : NOMS DES MAILLES SUPPORTS DES ELEMENTS
!     LOC   : LOCALISATION DES VALEURS (ELNO =>57, ELGA=>56,ELEM=>56)
!     CELD  : DESCRIPTEUR DU CHAM_ELEM (MODES LOCAUX,ADRESSES->.CELV)
!     NBNOMA: NOMBRE DE NOEUDS DE CHAQUE MAILLE
!     PERMUT: TABLEAU DES PERMUTATIONS DES NOEUDS DE CHAQUE TYPE-MA
!     MAXNOD: NOMBRE MAXI DE NOEUDS POUR LES DIFF. TYPE_MAILLES
!     TYPMA : TYPE_MAILLES
!     IR    : NUMERO D'ORDRE DU CHAMP
!     NBMAT : NOMBRE DE MAILLES A IMPRIMER
!     NUMMAI: NUMEROS DES MAILLES A IMPRIMER
!     LMASU : INDIQUE SI MAILLAGE ISSU DE SUPERTAB  .TRUE.
!     NBCMP : NOMBRE DE COMPOSANTES A IMPRIMER
!
!     ------------------------------------------------------------------
    character(len=3) :: toto
    character(len=4) :: nomgs
    character(len=8) :: nocmp, ktype
    character(len=24) :: nomst
    character(len=80) :: entete(10), titre, texte
    integer(kind=8) :: nbchs, nbcmpt, entier, nbspt, nnoe, ilong, imodel
    integer(kind=8) :: impre, iente, impel
    aster_logical :: afaire, lcmp, lnocen
!
!  --- INITIALISATIONS ----
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iachml, iad, iaec, iast, ibcmps, ic
    integer(kind=8) :: ichs, icm, icmax0, icmp, icmpl, icms
    integer(kind=8) :: icmsup, ico, icoef, icomax, icomm, icou, icp
    integer(kind=8) :: ida, idebu, idern, iel, ielg, ier, ies
    integer(kind=8) :: ifin, igre, igrel, ilig, imai
    integer(kind=8) :: inoa, inos
    integer(kind=8) :: ipg, ipoin1, ipoin2, ir, ires, irvg
    integer(kind=8) :: irvn, is0, isp, ispt, isup
    integer(kind=8) :: itseg2, itype, iutil, j, jj, jmax, jmod
    integer(kind=8) :: jspt, jt, jtitr, k, l, ll
    integer(kind=8) :: mode, nbcou, nbdats, nbelgr, nbpg, ncmpg
    integer(kind=8) :: ncmpp, nec, ni, npcalc, nsca, nscal
    integer(kind=8), pointer :: cmp_grel(:) => null()
    integer(kind=8), pointer :: ipcmps(:) => null()
    aster_logical, pointer :: ltabl(:) => null()
    integer(kind=8), pointer :: nbcmps_grel(:) => null()
    integer(kind=8), pointer :: nbcmpt_grel(:) => null()
    character(len=8), pointer :: nomchs(:) => null()
    character(len=8), pointer :: nomgds(:) => null()
    integer(kind=8), pointer :: perm(:) => null()
    integer(kind=8), pointer :: pos(:) => null()
    integer(kind=8), pointer :: scmp_dats(:) => null()
    integer(kind=8), pointer :: snbcps(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
    AS_ALLOCATE(vk8=nomgds, size=ncmpmx)
    AS_ALLOCATE(vk8=nomchs, size=ncmpmx)
    call wkvect('&&IRCERS.NBCMPS', 'V V I', ncmpmx, ibcmps)
    AS_ALLOCATE(vi=ipcmps, size=ncmpmx*ncmpmx)
    AS_ALLOCATE(vl=ltabl, size=ncmpmx)
    call jeveuo('&CATA.TE.MODELOC', 'L', imodel)
    call jeveuo(jexatr('&CATA.TE.MODELOC', 'LONCUM'), 'L', ilong)
!
    nomst = '&&IRECRI.SOUS_TITRE.TITR'
    call jeveuo(nomst, 'L', jtitr)
    titre = zk80(jtitr)
    do i = 1, ncmpmx
        ltabl(i) = .false.
    end do
    lnocen = .false.
    lcmp = .false.
!
!  --- RECHERCHE DES GRANDEURS SUPERTAB ----
!
    call irgags(ncmpmx, nomcmp, nomsym, nbchs, nomchs, &
                zi(ibcmps), nomgds, ipcmps)
!
!
!     NOM DE LA GRANDEUR SUPERTAB
    if (nbcmp .ne. 0) then
        lcmp = .true.
        if (nomgd .ne. 'VARI_R') then
            do i = 1, nbchs
                do j = 1, zi(ibcmps+i-1)
                    if (ncmps(1) .eq. ipcmps((i-1)*ncmpmx+j)) goto 899
                end do
            end do
899         continue
            nomgs = nomgds(i)
        end if
    end if
! --- DETERMINATION DU NOMBRE MAXIMUM DE SOUS-POINTS ---
    icomax = 0
    do igre = 1, nbgrel
        icoef = max(1, celd(4))
        if (icoef .gt. icomax) icomax = icoef
    end do
    icomm = 6
    if (ncmpu .eq. 0) then
        icmax0 = icomax
    else
        icmax0 = ncmpu
    end if
!
! -- ALLOCATION DES TABLEAUX DE TRAVAIL ---
!
    if (loc .eq. 'ELNO') then
        call jedetr('&&IRCERS.VALNOE')
        call wkvect('&&IRCERS.VALNOE', 'V V R', ncmpmx*icomm, irvn)
    else if (loc .eq. 'ELGA' .or. loc .eq. 'ELEM') then
        call jedetr('&&IRCERS.VALGAU')
        call wkvect('&&IRCERS.VALGAU', 'V V R', ncmpmx*icomm, irvg)
    end if
!
! ---- BOUCLE SUR LES DIVERSES GRANDEURS SUPERTAB ----
    if (nbcmp .ne. 0 .and. nomgd .ne. 'VARI_R') nbchs = 1
    do ichs = 1, nbchs
        if (ichs .gt. 1) then
            afaire = .false.
            do icp = 1, zi(ibcmps-1+ichs)
                afaire = (afaire .or. ltabl((ipcmps((ichs-1)* &
                                                    ncmpmx+icp))))
            end do
            if (.not. afaire) goto 10
        end if
!
        call jedetr('&&IRCERS.SOUS_PT')
        call wkvect('&&IRCERS.SOUS_PT', 'V V I', icomax, jspt)
!
! --- GROUPEMENT DES VARIABLES SCALAIRES 6 PAR 6 ----
!  ---  DETERMINATION DU NOMBRE DE DATASETS SUPERTAB A IMPRIMER --
!
        if (zi(ibcmps-1+ichs) .eq. 1 .and. icomax .gt. 1) then
            do i = 1, icmax0
                zi(jspt-1+i) = 6
            end do
            ilig = icmax0/6
            ires = icmax0-ilig*6
            if (ires .eq. 0) then
                nbdats = ilig
            else
                nbdats = ilig+1
                zi(jspt-1+nbdats) = ires
            end if
            nbcmpt = 6
            nomgs = 'VARI'
        else if (nbcmp .ne. 0 .and. nomgd .ne. 'VARI_R') then
            AS_ALLOCATE(vi=scmp_dats, size=nbcmp)
            ilig = nbcmp/6
            ires = nbcmp-ilig*6
            ni = 0
            scmp_dats(1) = ni
            if (ires .eq. 0) then
                nbdats = ilig
                do i = 1, nbdats
                    zi(ibcmps+i-1) = 6
                    ni = ni+6
                    scmp_dats(1+i) = ni
                end do
            else
                nbdats = ilig+1
                do i = 1, nbdats-1
                    zi(ibcmps+i-1) = 6
                    ni = ni+6
                    scmp_dats(1+i) = ni
                end do
                zi(ibcmps+nbdats-1) = ires
                scmp_dats(nbdats+1) = ni+ires
            end if
            nbcmpt = 6
        else
            nbdats = icomax
            do i = 1, nbdats
                zi(jspt-1+i) = 1
            end do
            nbcmpt = zi(ibcmps-1+ichs)
            nomgs = nomgds(ichs)
        end if
!
!
!
! --- ECRITURE DE L'ENTETE SUPERTAB ----
!
        call ecrtes(nomsd, titr, nomgs, ir, loc, &
                    nbcmpt, 2, entete, lcmp)
!
!
! --- POSITION DES COMPOSANTES SELECTIONNEES ---
!
        if (nbcmp .ne. 0 .and. nomgd .ne. 'VARI_R') then
            k = 0
            l = 0
            ll = 0
            AS_ALLOCATE(vi=perm, size=nbcmp*nbgrel)
            AS_ALLOCATE(vi=nbcmps_grel, size=nbgrel)
            AS_ALLOCATE(vi=pos, size=nbcmp*nbgrel)
            AS_ALLOCATE(vi=nbcmpt_grel, size=nbgrel)
            AS_ALLOCATE(vi=snbcps, size=nbgrel+1)
            AS_ALLOCATE(vi=cmp_grel, size=nbcmp*nbgrel)
            snbcps(1) = 0
            do igrel = 1, nbgrel
                mode = celd(celd(4+igrel)+2)
                ipoin1 = longr(igrel)
                ipoin2 = longr(igrel+1)
                nbelgr = ipoin2-ipoin1-1
                if (mode .eq. 0) goto 904
                jmod = imodel+zi(ilong-1+mode)-1
                nec = nbec(zi(jmod-1+2))
                call jedetr('&&IRCERS.ENT_COD')
                call wkvect('&&IRCERS.ENT_COD', 'V V I', nec, iaec)
                call dgmode(mode, imodel, ilong, nec, zi(iaec))
                ncmpg = 0
!             POSITIONS DES COMPOSANTES SELECTIONNEES PRESENTES
!             DANS LE  GREL PARMI LES COMPOSANTES SELECTIONNEES
                do icmpl = 1, nbcmp
                    if (exisdg(zi(iaec), ncmps(icmpl))) then
                        cmp_grel(1+k) = ncmps(icmpl)
                        k = k+1
                        ncmpg = ncmpg+1
                    end if
                end do
!             SOMME DES COMPOSANTES SELECTIONNEES PAR GREL
                snbcps(igrel+1) = k
!             NOMBRE DE COMPOSANTES SELECTIONNEES PRESENTES
!             DANS LE GREL
                nbcmps_grel(igrel) = ncmpg
                if (igrel .eq. nbgrel .and. k .lt. nbcmp) then
                    call utmess('F', 'PREPOST_83')
                end if
                ncmpp = 0
                do i = 1, ncmpmx
                    if (exisdg(zi(iaec), i)) then
                        ncmpp = ncmpp+1
!                   POSITIONS DES COMPOSANTES SELECTIONNEES PRESENTES
!                   DANS LE GREL PARMI LES COMPOSANTES DU GREL
                        do j = 1, nbcmps_grel(igrel)
                            if (i .eq. cmp_grel(1+j-1+snbcps(igrel))) then
                                pos(1+l) = ncmpp
                                l = l+1
                            end if
                        end do
!                   POSITION DES COMPOSANTES SELECTIONNEES
!                   DANS LES DATASETS
                        do j = 1, nbcmp
                            if (i .eq. ncmps(j)) then
                                perm(ll+1) = j
                                ll = ll+1
                            end if
                        end do
                    end if
                end do
!             NOMBRE DE COMPOSANTES DANS LE GREL
                nbcmpt_grel(igrel) = ncmpp
904             continue
            end do
        end if
!
!
! --- IMPRESSION DES DATASETS SUPERTAB ---
!
!
        do ida = 1, nbdats
!
            iente = 1
            impre = 0
            ifin = 1
            idebu = 1
            entete(4) = ' '
            texte = ' '
!
! --- ECRITURE DANS L'ENTETE SUPERTAB DES NOMS DE COMPOSANTES---
!
            if (nbcmp .ne. 0 .and. nomgd .ne. 'VARI_R') then
                do icp = 1, zi(ibcmps+ida-1)
                    nocmp = nocmpl(scmp_dats(ida)+icp)
                    iutil = lxlgut(nocmp)
                    ifin = idebu+iutil
                    texte(idebu:ifin) = nocmp(1:iutil)//' '
                    idebu = ifin+1
                end do
                texte(ifin+2:ifin+7) = '('//loc//')'
            else
                nbspt = zi(jspt-1+ida)
                do icp = 1, zi(ibcmps-1+ichs)
                    if (zi(ibcmps-1+ichs) .eq. 1 .and. icomax .gt. 1) then
                        do ispt = 1, nbspt
                            if (ncmpu .eq. 0) then
                                entier = (ida-1)*6+ispt
                            else
                                entier = nucmp((ida-1)*6+ispt)
                            end if
                            nocmp = nomcmp(ipcmps((ichs-1)*ncmpmx+icp))
                            iutil = lxlgut(nocmp)
                            call codent(entier, 'G', toto)
                            ifin = idebu+iutil+3
                            texte(idebu:ifin) = nocmp(1:iutil)//'_'// &
                                                toto//' '
                            idebu = ifin+1
                        end do
                    else
                        nocmp = nomcmp(ipcmps((ichs-1)*ncmpmx+icp))
                        iutil = lxlgut(nocmp)
                        ifin = idebu+iutil
                        texte(idebu:ifin) = nocmp(1:iutil)//' '
                        idebu = ifin+1
                    end if
                end do
                texte(ifin+2:ifin+7) = '('//loc
                idern = ifin+7
                if (zi(ibcmps-1+ichs) .gt. 1 .and. icomax .gt. 1) then
                    call codent(ida, 'G', toto)
                    texte(ifin+8:ifin+14) = 'SPT_'//toto
                    idern = ifin+14
                end if
                texte(idern:idern+1) = ')'
            end if
!
            iutil = lxlgut(texte)
            jmax = lxlgut(titre)
            jmax = min(jmax, (80-iutil-3))
            entete(4) = titre(1:jmax)//' - '//texte(1:iutil)
!
            do igrel = 1, nbgrel
                mode = celd(celd(4+igrel)+2)
                ipoin1 = longr(igrel)
                ipoin2 = longr(igrel+1)
                nbelgr = ipoin2-ipoin1-1
                if (mode .eq. 0) goto 12
                jmod = imodel+zi(ilong-1+mode)-1
                nec = nbec(zi(jmod-1+2))
                call jedetr('&&IRCERS.ENT_COD')
                call wkvect('&&IRCERS.ENT_COD', 'V V I', nec, iaec)
                call dgmode(mode, imodel, ilong, nec, zi(iaec))
                iad = celd(celd(4+igrel)+8)
                nscal = digdel(mode)
!
                if (nbcmp .ne. 0 .and. nomgd .ne. 'VARI_R') then
                    nsca = nscal
                    if (nbcmps_grel(igrel) .ne. 0) then
                        goto 62
                    else
                        goto 12
                    end if
                end if
!
                icoef = max(1, celd(4))
                if (zi(ibcmps-1+ichs) .eq. 1 .and. icomax .gt. 1) then
                    if (ncmpu .eq. 0) then
                        ico = (ida-1)*6+1
                    else
                        ico = 1
                    end if
                else
                    ico = ida
                end if
                if (icoef .lt. ico) goto 12
                nsca = nscal*icoef
                ncmpp = 0
                do i = 1, ncmpmx
                    if (exisdg(zi(iaec), i)) then
                        ncmpp = ncmpp+1
                        if (ichs .eq. 1) ltabl(i) = .true.
                    end if
                end do
                do i = 1, zi(ibcmps-1+ichs)
                    if (exisdg(zi(iaec), ipcmps((ichs-1)*ncmpmx+i))) goto 62
                end do
                goto 12
62              continue
                do ielg = 1, nbelgr
                    iel = ligrel(ipoin1+ielg-1)
                    if (iel .le. 0) goto 13
!
! --- IMPRESSION DU CHAMELEM SUR UNE LISTE DE MAILLES ---
!
                    if (nbmat .ne. 0) then
                        do imai = 1, nbmat
                            if (iel .eq. nummai(imai)) goto 15
                        end do
                        goto 13
                    end if
15                  continue
                    impel = 1
!
!           RECHERCHE DE L'ADRESSE DANS VALE DU DEBUT DES VALEURS
!
                    iachml = iad+nsca*(ielg-1)
!
!
!    --- CHAMELEM AUX NOEUDS ET AU POINTS DE GAUSS
!
!       - ELEMENTS NON DISPONIBLES DANS IDEAS
!
                    itype = typma(iel)
                    call jenuno(jexnum('&CATA.TM.NOMTM', itype), ktype)
                    if (ktype .eq. 'PYRAM5') then
                        call utmess('A', 'PREPOST_84')
                    else if (ktype .eq. 'PYRAM13') then
                        call utmess('A', 'PREPOST_85')
                    else
!
!    --- CHAMELEM AUX NOEUDS ---
!
                        if (loc .eq. 'ELNO') then
!
                            if (nbcmp .ne. 0 .and. nomgd .ne. 'VARI_R') then
!
                                npcalc = nscal/nbcmpt_grel(igrel)
                                nnoe = nbnoma(iel)
                                itype = typma(iel)
                                call jenuno(jexnum('&CATA.TM.NOMTM', itype), ktype)
                                if (ktype .eq. 'HEXA27') then
                                    nnoe = nnoe-7
                                    lnocen = .true.
                                else if (ktype .eq. 'TRIA7') then
                                    nnoe = nnoe-1
                                    lnocen = .true.
                                else if (ktype .eq. 'TETRA15') then
                                    nnoe = nnoe-5
                                    lnocen = .true.
                                else if (ktype .eq. 'PENTA18') then
                                    nnoe = nnoe-3
                                    lnocen = .true.
                                else if (ktype .eq. 'PENTA21') then
                                    nnoe = nnoe-6
                                    lnocen = .true.
                                else if (ktype .eq. 'QUAD9') then
                                    nnoe = nnoe-1
                                    lnocen = .true.
                                else if (ktype .eq. 'SEG4') then
                                    nnoe = nnoe-2
                                    call jenonu(jexnom('&CATA.TM.NOMTM', 'SEG2'), itseg2)
                                    itype = itseg2
                                end if
                                nbcou = npcalc/nnoe
!
                                do inos = 1, nnoe
                                    inoa = 0
                                    do iast = 1, nnoe
                                        isup = permut(iast, itype)
                                        if (inos .eq. isup) then
                                            inoa = iast
                                            goto 529
                                        end if
                                    end do
!
529                                 continue
                                    ASSERT(inoa .ne. 0)
!
                                    do icou = 1, nbcou
                                        jj = iachml-1+nbcmpt_grel(igrel) &
                                             *(inoa-1)+(icou-1)*nbcmpt_grel(1+ &
                                                                            igrel-1)*nnoe
!
                                        do i = 1, nbcmpt
                                            zr(irvn-1+i) = 0.d0
                                        end do
!
                                        do icm = 1, nbcmps_grel(igrel)
                                            j = perm(1+snbcps(igrel)+ &
                                                     icm-1)
                                            if (j .le. 6*ida .and. j .ge. (6*ida-5)) then
                                                jt = j-6*(ida-1)
                                                ic = pos(1+snbcps(igrel)+ &
                                                         icm-1)
                                                zr(irvn-1+jt) = vale(jj+ic)
                                            end if
                                        end do
!
                                        if (iente .eq. 1) then
                                            write (ifi, '(A80)') (entete(i), &
                                                                  i=1, 10)
                                            iente = 0
                                        end if
                                        if (impel .eq. 1) then
                                            if (lmasu) then
                                                call lxliis(nomel(iel) (2:8), ies, ier)
                                            else
                                                ies = iel
                                            end if
                                            write (ifi, '(4I10,5X,A,A)') &
                                                ies, 1, nnoe, nbcmpt, '% MAILLE ', &
                                                nomel(iel)
                                            impel = 0
                                        end if
                                        write (ifi, '(6(1PE13.5E3))') ( &
                                            zr(irvn-1+i), i=1, nbcmpt)
!
                                    end do
                                end do
!
                            else
!
                                npcalc = nscal/ncmpp
                                nnoe = nbnoma(iel)
                                itype = typma(iel)
                                call jenuno(jexnum('&CATA.TM.NOMTM', itype), ktype)
                                if (ktype .eq. 'HEXA27') then
                                    nnoe = nnoe-7
                                    lnocen = .true.
                                else if (ktype .eq. 'TRIA7') then
                                    nnoe = nnoe-1
                                    lnocen = .true.
                                else if (ktype .eq. 'PENTA18') then
                                    nnoe = nnoe-3
                                    lnocen = .true.
                                else if (ktype .eq. 'QUAD9') then
                                    nnoe = nnoe-1
                                    lnocen = .true.
                                else if (ktype .eq. 'SEG4') then
                                    nnoe = nnoe-2
                                    call jenonu(jexnom('&CATA.TM.NOMTM', 'SEG2'), itseg2)
                                    itype = itseg2
                                end if
                                nbcou = npcalc/nnoe
!
                                do inos = 1, nnoe
                                    inoa = 0
                                    do iast = 1, nnoe
                                        isup = permut(iast, itype)
                                        if (inos .eq. isup) then
                                            inoa = iast
                                            goto 29
                                        end if
                                    end do
29                                  continue
                                    ASSERT(inoa .ne. 0)
                                    do icou = 1, nbcou
                                        j = iachml-1+ncmpp*icoef*(inoa- &
                                                                  1)+(icou-1)*ncmpp*icoef*nnoe+ &
                                            ncmpp*(ico-1)
                                        do i = 1, nbcmpt
                                            zr(irvn-1+i) = 0.d0
                                        end do
                                        ic = 0
                                        do icmp = 1, ncmpmx
                                            if (exisdg(zi(iaec), icmp)) then
                                                ic = ic+1
                                                do icms = 1, zi(ibcmps-1+ichs)
                                                    icmsup = ipcmps((ichs-1)*ncmpmx+icms)
                                                    if (icmp .eq. icmsup) then
                                                        impre = 1
                                                        do isp = 1, zi(jspt-1+ida)
                                                            zr(irvn-1+icms-1+isp) = &
                                                                vale(j+ic+ncmpp*(isp-1))
                                                        end do
                                                        goto 22
                                                    end if
                                                end do
                                            end if
22                                          continue
                                        end do
                                        if (impre .eq. 1) then
                                            if (iente .eq. 1) then
                                                write (ifi, '(A80)') (entete(i), &
                                                                      i=1, 10)
                                                iente = 0
                                            end if
                                            if (impel .eq. 1) then
                                                if (lmasu) then
                                                    call lxliis(nomel(iel) (2:8), ies, ier)
                                                else
                                                    ies = iel
                                                end if
                                                write (ifi, '(4I10,5X,A,A)') &
                                                    ies, 1, nnoe, nbcmpt, '% MAILLE ', &
                                                    nomel(iel)
                                                impel = 0
                                            end if
                                            write (ifi, '(6(1PE13.5E3))') ( &
                                                zr(irvn-1+i), i=1, nbcmpt)
                                        end if
                                    end do
                                end do
                            end if
!
!  --- CHAMELEM AUX POINTS DE GAUSS---
!
                        else if (loc .eq. 'ELGA' .or. loc .eq. 'ELEM') then
!
                            if (nbcmp .ne. 0 .and. nomgd .ne. 'VARI_R') then
                                nbpg = nscal/nbcmpt_grel(igrel)
!
                                do i = 1, nbcmpt
                                    zr(irvg-1+i) = 0.d0
                                end do
!
                                do icm = 1, nbcmps_grel(igrel)
                                    j = perm(1+snbcps(igrel)+icm- &
                                             1)
!
                                    if (j .le. 6*ida .and. j .ge. (6*ida-5)) then
                                        jt = j-6*(ida-1)
!
                                        do ipg = 1, nbpg
                                            jj = nbcmpt_grel(igrel)*(ipg-1)+ &
                                                 pos(1+snbcps(igrel)+ &
                                                     icm-1)
                                            zr(irvg-1+jt) = zr(irvg-1+jt)+ &
                                                            vale(iachml-1+jj)
                                        end do
                                        zr(irvg-1+jt) = zr(irvg-1+jt)/ &
                                                        nbpg
!
                                    end if
                                end do
!
                                if (iente .eq. 1) then
                                    write (ifi, '(A80)') (entete(i), i=1, &
                                                          10)
                                    iente = 0
                                end if
                                if (lmasu) then
                                    call lxliis(nomel(iel) (2:8), ies, ier)
                                else
                                    ies = iel
                                end if
                                write (ifi, '(2I10,5X,2A)') ies, nbcmpt, &
                                    '% MAILLE ', nomel(iel)
                                write (ifi, '(6(1PE13.5E3))') (zr(irvg- &
                                                                  1+i), i=1, nbcmpt)
!
!
                            else
                                npcalc = nscal/ncmpp
                                nbpg = npcalc
                                do i = 1, nbcmpt
                                    zr(irvg-1+i) = 0.d0
                                end do
                                ic = 0
                                do icmp = 1, ncmpmx
                                    if (exisdg(zi(iaec), icmp)) then
                                        ic = ic+1
                                        do icms = 1, zi(ibcmps-1+ichs)
                                            icmsup = ipcmps((ichs-1)*ncmpmx+icms)
                                            if (icmp .eq. icmsup) then
                                                impre = 1
                                                do isp = 1, zi(jspt-1+ida)
                                                    if (ncmpu .eq. 0) then
                                                        is0 = isp
                                                    else
                                                        is0 = nucmp((ida-1)*6+isp)
                                                    end if
                                                    do ipg = 1, nbpg
!
                                                        j = iachml-1+ncmpp*icoef*(ipg-1) &
                                                            +ncmpp*(ico-1)
                                                        zr(irvg-1+icms-1+isp) = &
                                                            zr(irvg-1+icms-1+isp)+ &
                                                            vale(j+ic+ncmpp*(is0-1))
!
                                                    end do
                                                    zr(irvg-1+icms-1+isp) = zr(irvg- &
                                                                               1+icms-1+isp)/nbpg
!
                                                end do
                                                goto 19
                                            end if
                                        end do
                                    end if
19                                  continue
                                end do
                                if (impre .eq. 1) then
                                    if (iente .eq. 1) then
                                        write (ifi, '(A80)') (entete(i), &
                                                              i=1, 10)
                                        iente = 0
                                    end if
                                    if (impel .eq. 1) then
                                        if (lmasu) then
                                            call lxliis(nomel(iel) (1:8), ies, ier)
                                        else
                                            ies = iel
                                        end if
                                        write (ifi, '(2I10,5X,2A)') ies, &
                                            nbcmpt, '% MAILLE ', nomel(iel)
                                        impel = 0
                                    end if
                                    write (ifi, '(6(1PE13.5E3))') (zr( &
                                                                   irvg-1+i), i=1, nbcmpt)
                                    impre = 0
                                end if
                            end if
                        end if
                    end if
13                  continue
                end do
12              continue
            end do
            if (iente .eq. 0) write (ifi, '(A)') '    -1'
        end do
10      continue
    end do
!
    if (lnocen) then
        call utmess('A', 'PREPOST_86')
    end if
!
    call jedetr('&&IRCERS.VALNOE')
    call jedetr('&&IRCERS.VALGAU')
    call jedetr('&&IRCERS.SOUS_PT')
    call jedetr('&&IRCERS.ENT_COD')
    AS_DEALLOCATE(vk8=nomgds)
    AS_DEALLOCATE(vk8=nomchs)
    call jedetr('&&IRCERS.NBCMPS')
    AS_DEALLOCATE(vi=ipcmps)
    AS_DEALLOCATE(vl=ltabl)
    AS_DEALLOCATE(vi=perm)
    AS_DEALLOCATE(vi=nbcmps_grel)
    AS_DEALLOCATE(vi=pos)
    AS_DEALLOCATE(vi=nbcmpt_grel)
    AS_DEALLOCATE(vi=snbcps)
    AS_DEALLOCATE(vi=cmp_grel)
    AS_DEALLOCATE(vi=scmp_dats)
    call jedema()
end subroutine
