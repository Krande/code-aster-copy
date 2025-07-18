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
! person_in_charge: nicolas.sellenet at edf.fr
! aslint: disable=W1504
!
subroutine lrcame(nrofic, nochmd, nomamd, nomaas, ligrel, &
                  option, param, typech, typen, npgma, &
                  npgmm, nspmm, nbcmpv, ncmpva, ncmpvm, &
                  iinst, numpt, numord, inst, crit, &
                  prec, nomgd, ncmprf, jnocmp, chames, &
                  codret)
!
    use as_med_module, only: as_med_open
    implicit none
!
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_mficlo.h"
#include "asterfort/as_mficom.h"
#include "asterfort/as_mfinvr.h"
#include "asterfort/as_mlbnuv.h"
#include "asterfort/as_mmhnme.h"
#include "asterfort/assert.h"
#include "asterfort/cescre.h"
#include "asterfort/cnscre.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lrcmle.h"
#include "asterfort/lrcmpr.h"
#include "asterfort/lrcmva.h"
#include "asterfort/lrcmve.h"
#include "asterfort/lrmpga.h"
#include "asterfort/lrmtyp.h"
#include "asterfort/mdchin.h"
#include "asterfort/mdexch.h"
#include "asterfort/mdexma.h"
#include "asterfort/mdexpm.h"
#include "asterfort/ulisog.h"
#include "asterfort/utlicm.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nrofic, typen
    integer(kind=8) :: ncmprf, jnocmp
    integer(kind=8) :: nbcmpv
    integer(kind=8) :: iinst, numpt, numord
    integer(kind=8) :: npgma(*), npgmm(*), nspmm(*)
    integer(kind=8) :: codret, codre2
    character(len=*) :: typech
    character(len=8) :: nomgd, nomaas
    character(len=8) :: crit, param
    character(len=19) :: chames, ligrel
    character(len=24) :: option
    character(len=*) :: nochmd, nomamd
    character(len=*) :: ncmpva, ncmpvm
    real(kind=8) :: inst, prec
!
! --------------------------------------------------------------------------------------------------
!
!     LECTURE D'UN CHAMP - FORMAT MED
!
! --------------------------------------------------------------------------------------------------
!
!      ENTREES:
!        NROFIC : UNITE LOGIQUE DU FICHIER MED
!        NOCHMD : NOM MED DU CHAMP A LIRE
!        NOMAMD : NOM MED DU MAILLAGE LIE AU CHAMP A LIRE
!                  SI ' ' : ON SUPPOSE QUE C'EST LE PREMIER MAILLAGE
!                          DU FICHIER
!        NOMAAS : NOM ASTER DU MAILLAGE
!        NBVATO : NOMBRE DE VALEURS TOTAL
!        TYPECH : TYPE DE CHAMP AUX ELEMENTS : ELEM/ELGA/ELNO/NOEU
!        TYPEN  : TYPE D'ENTITE DU CHAMP
!                (MED_NOEUD=3,MED_MAILLE=0,MED_NOEUD_MAILLE=4)
!        NPGMA  : NOMBRE DE POINTS DE GAUSS PAR MAILLE (ASTER)
!        NPGMM  : NOMBRE DE POINTS DE GAUSS PAR MAILLE (MED)
!        NSPMM  : NOMBRE DE SOUS-POINTS PAR MAILLE (MED)
!        NBCMPV : NOMBRE DE COMPOSANTES VOULUES
!                 SI NUL, ON LIT LES COMPOSANTES A NOM IDENTIQUE
!        NCMPVA : LISTE DES COMPOSANTES VOULUES POUR ASTER
!        NCMPVM : LISTE DES COMPOSANTES VOULUES DANS MED
!        IINST  : 1 SI LA DEMANDE EST FAITE SUR UN INSTANT, 0 SINON
!        NUMPT  : NUMERO DE PAS DE TEMPS EVENTUEL
!        NUMORD : NUMERO D'ORDRE EVENTUEL DU CHAMP
!        INST   : INSTANT EVENTUEL
!        CRIT   : CRITERE SUR LA RECHERCHE DU BON INSTANT
!        PREC   : PRECISION SUR LA RECHERCHE DU BON INSTANT
!        NOMGD  : NOM DE LA GRANDEUR ASSOCIEE AU CHAMP
!        NCMPRF : NOMBRE DE COMPOSANTES DE REFERENCE DU CHAMP SIMPLE
!        JNOCMP : ADRESSE DU NOM DES COMP. DE REF. DU CHAMP SIMPLE
!      SORTIES:
!         CHAMES : NOM DU CHAMP A CREER
!         CODRET : CODE DE RETOUR (0 : PAS DE PB, NON NUL SI PB)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=6), parameter :: nompro = 'LRCAME'
    integer(kind=8) :: typent, vali(4)
    character(len=64), parameter :: ednopf = ' '
    integer(kind=8), parameter :: ednoeu = 3, edmail = 0, edconn = 1, ednoda = 0
    integer(kind=8), parameter :: edlect = 0, typnoe = 0, ntypel = 26, npgmax = 27
    integer(kind=8) :: iaux, letype, vlib(3), vfic(3), iret
    med_idt :: idfimd, ifimed
    integer(kind=8) :: indpg(ntypel, npgmax)
    integer(kind=8) :: nbvato, nbcmfi, nbval, nbma, nbprof
    integer(kind=8) :: adsl, adsv, adsd, i, j
    integer(kind=8) :: ncmput, existc, ndim, npas, lgproa
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: tygeom, nbtyp
    integer(kind=8) :: nnotyp(MT_NTYMAX), modnum(MT_NTYMAX), numnoa(MT_NTYMAX, MT_NNOMAX)
    integer(kind=8) :: typgeo(MT_NTYMAX), lygeom(MT_NTYMAX), lypent(MT_NTYMAX), ltyp(MT_NTYMAX)
    integer(kind=8) :: renumd(MT_NTYMAX), nlyval(MT_NTYMAX), nuanom(MT_NTYMAX, MT_NNOMAX)
    integer(kind=8) :: nbtylu, iaux2, k, nbty(MT_NTYMAX), lnbpro(MT_NTYMAX)
    integer(kind=8) :: nbnoma, nmatyp, jntpro, lgprof, cptyma, iprof
    integer(kind=8) :: jnumty, numma, ima, hdfok, medok, jmaill
    aster_logical :: lrenum
    character(len=1) :: saux01
    character(len=8) :: saux08, modele
    character(len=8) :: nomtyp(MT_NTYMAX)
    character(len=19) :: prefix
    character(len=24) :: numcmp, ntncmp, ntucmp, ntvale, nmcmfi(MT_NTYMAX)
    character(len=64) :: valk(2)
    character(len=24) :: ntproa, nmcmfl
    character(len=64) :: nomprf
    character(len=200) :: nofimd
    character(len=255) :: kfic
    character(len=2) :: k2bid
    real(kind=8) :: valr
    aster_logical :: existm, existt
    aster_logical :: logaux
    character(len=8), pointer :: lgrf(:) => null()
    integer(kind=8), pointer :: nume(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
    real(kind=8), pointer :: vinst(:) => null()
    integer(kind=8) :: iCmp
    character(len=8), pointer :: cmpUserName(:) => null()
    character(len=8), pointer :: cmpCataName(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call infniv(ifm, niv)
    nomprf = ' '
!
    if (niv .gt. 1) then
        write (ifm, 101) 'DEBUT DE '//nompro
        write (ifm, *) '.. NOM DU CHAMP A LIRE : ', nochmd
    end if
101 format(/, 10('='), a, 10('='),/)
!
! 1.2. ==> NOMS DES TABLEAUX DE TRAVAIL
!
    numcmp = '&&'//nompro//'.NUMERO_CMP     '
    ntncmp = '&&'//nompro//'.NOMCMP         '
    ntucmp = '&&'//nompro//'.UNITECMP       '
    ntvale = '&&'//nompro//'.VALEUR         '
    ntproa = '&&'//nompro//'.PROFIL_ASTER   '
    prefix = '&&'//nompro//'.MED'
!
! 1.3. ==> NOM DU FICHIER MED
!
    call ulisog(nrofic, kfic, saux01)
    if (kfic(1:1) .eq. ' ') then
        call codent(nrofic, 'G', saux08)
        nofimd = 'fort.'//saux08
    else
        nofimd = kfic(1:200)
    end if
!
    if (niv .gt. 1) then
        write (ifm, *) '<', nompro, '> NOM DU FICHIER MED : ', nofimd
    end if
!
! 1.4. ==> VERIFICATION DU FICHIER MED
!
! 1.4.1. ==> VERIFICATION DE LA VERSION HDF
!
    call as_mficom(nofimd, hdfok, medok, codret)
    if (hdfok .eq. 0) then
        saux08 = 'mficom'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
!
! 1.4.2. ==> VERIFICATION DE LA VERSION MED
!
    if (medok .eq. 0) then
        vali(1) = codret
        call utmess('F+', 'MED_24')
        call as_mlbnuv(vlib(1), vlib(2), vlib(3), iret)
        if (iret .eq. 0) then
            vali(1) = vlib(1)
            vali(2) = vlib(2)
            vali(3) = vlib(3)
            call utmess('F+', 'MED_25', ni=3, vali=vali)
        end if
        call as_med_open(idfimd, nofimd, edlect, codret)
        call as_mfinvr(idfimd, vfic(1), vfic(2), vfic(3), iret)
        if (iret .eq. 0) then
            if (vfic(2) .eq. -1 .or. vfic(3) .eq. -1) then
                call utmess('F+', 'MED_26')
            else
                vali(1) = vfic(1)
                vali(2) = vfic(2)
                vali(3) = vfic(3)
                call utmess('F+', 'MED_27', ni=3, vali=vali)
            end if
            if (vfic(1) .lt. vlib(1) .or. (vfic(1) .eq. vlib(1) .and. vfic(2) .lt. vlib(2)) &
                .or. &
                ( &
                vfic(1) .eq. vlib(1) .and. vfic(2) .eq. vlib(2) .and. vfic(3) .eq. vlib(3) &
                )) then
                call utmess('F+', 'MED_28')
            end if
        end if
        call as_mficlo(idfimd, codret)
    end if
!
! 1.5. ==> VERIFICATION DE L'EXISTENCE DU MAILLAGE CONCERNE
!
! 1.5.1. ==> C'EST LE PREMIER MAILLAGE DU FICHIER
!            ON RECUPERE SON NOM ET SA DIMENSION.
!
    if (nomamd .eq. ' ') then
        ifimed = 0
        call mdexpm(nofimd, ifimed, nomamd, existm, ndim, &
                    codret)
        if (.not. existm) then
            call utmess('F', 'MED_50', sk=nofimd)
        end if
!
! 1.5.2. ==> C'EST UN MAILLAGE DESIGNE PAR UN NOM
!            ON RECUPERE SA DIMENSION.
!
    else
        iaux = 1
        ifimed = 0
        call mdexma(nofimd, ifimed, nomamd, iaux, existm, &
                    ndim, codret)
        if (.not. existm) then
            valk(1) = nomamd(1:24)
            valk(2) = nofimd(1:24)
            call utmess('F', 'MED_51', nk=2, valk=valk)
        end if
    end if
!
    if (typech .eq. 'NOEU') then
        call dismoi('NB_NO_MAILLA', nomaas, 'MAILLAGE', repi=nbnoma)
        nbvato = nbnoma
    else
        call dismoi('NB_MA_MAILLA', nomaas, 'MAILLAGE', repi=nbma)
        nbvato = nbma
    end if
!
    if (niv .gt. 1) then
        write (ifm, *) '.. NOM DU MAILLAGE MED ASSOCIE : ', nomamd
        write (ifm, *) '   DE DIMENSION ', ndim
    end if
!
! 2.2. ==> VERIFICATIONS DES COMPOSANTES ASTER DEMANDEES
!          EN SORTIE, ON A :
!       NCMPUT : NOMBRE DE COMPOSANTES VALIDES.
!       NUMCMP : TABLEAU DES NUMEROS DES COMPOSANTES VALIDES
!       NTNCMP : TABLEAU DES NOMS DES COMPOSANTES VALIDES (K8)
!       NTUCMP : TABLEAU DES UNITES DES COMPOSANTES VALIDES (K16)
!
    if (nbcmpv .ne. 0) then
        call jeveuo(ncmpva, 'L', iaux)
    else
        iaux = 1
    end if
!
    if (nbcmpv > 0) then
        AS_ALLOCATE(vk8=cmpUserName, size=nbcmpv)
        do iCmp = 1, nbcmpv
            cmpUserName(iCmp) = zk8(iaux-1+iCmp)
        end do
    end if
    AS_ALLOCATE(vk8=cmpCataName, size=ncmprf)
    do iCmp = 1, ncmprf
        cmpCataName(iCmp) = zk8(jnocmp-1+iCmp)
    end do
    call utlicm(nomgd, &
                nbcmpv, cmpUserName, &
                ncmprf, cmpCataName, &
                ncmput, numcmp, &
                ntncmp, ntucmp)
    AS_DEALLOCATE(vk8=cmpUserName)
    AS_DEALLOCATE(vk8=cmpCataName)
!
!====
! 2. OUVERTURE DU FICHIER EN LECTURE
!====
!
    call as_med_open(idfimd, nofimd, edlect, codret)
    if (codret .ne. 0) then
        saux08 = 'mfiope'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
!
! 2.1. ==> . RECUPERATION DES NB/NOMS/NBNO/NBITEM DES TYPES DE MAILLES
!            DANS CATALOGUE
!          . RECUPERATION DES TYPES GEOMETRIE CORRESPONDANT POUR MED
!          . VERIF COHERENCE AVEC LE CATALOGUE
!
    call lrmtyp(nbtyp, nomtyp, nnotyp, typgeo, renumd, &
                modnum, nuanom, numnoa)
!
! 2.1.1 ==> LE CHAMP EXISTE-T-IL DANS LE FICHIER ?
!          AU BON NUMERO D'ORDRE ?
!          CONTIENT-IL A MINIMA LES COMPOSANTES VOULUES ?
!          LE NOMBRE DE VALEURS EST-IL CORRECT ?
!          SI OUI A TOUTES SES QUESTIONS, EXISTC VAUT 3.
!          ON RECUPERE ALORS :
!          . LE NOMBRE DE COMPOSANTES QU'IL Y A.
!          . LE NOM DE CES COMPOSANTES.
!
!     COMME DANS LRMMMA :
!    REMARQUE : GRACE A LA RENUMEROTATION, ON PARCOURT LES TYPES DE
!    MAILLES DANS L'ORDRE CROISSANT DE LEUR TYPE MED. CE N'EST PAS
!    OBLIGATOIRE SI ON DISPOSE DES TABLEAUX DE NUMEROTATION DES MAILLES.
!    MAIS QUAND CES TABLEAUX SONT ABSENTS, LES CONVENTIONS MED PRECISENT
!    QUE LA NUMEROTATION IMPLICITE SE FAIT DANS CET ORDRE. DONC ON
!    LE FAIT !
!
    nbtylu = 0
    nbcmfi = 0
    existt = ASTER_FALSE
!
    numma = 1
!     EN SORTIE DE MDEXMA/MDEXPM, CODRET=0
    codret = 0
    do letype = 0, nbtyp
!
! 2.2.1. ==> LES BONS TYPES
!
        iaux = letype
!
        if (iaux .eq. 0) then
            typent = ednoeu
            tygeom = typnoe
        else
            iaux = renumd(iaux)
            typent = typen
            tygeom = typgeo(iaux)
        end if
!
!       RECUPERE LE NOMBRE DE MAILLES DE TYPE TYGEOM

        ! incompatibilite (tygeom==0 <=> MED_NONE) et (edmail==0 <=> MED_CELL)
        ! => on saute pour eviter une Erreur d'appel de l'API dans MED
        if (tygeom .eq. 0 .and. edmail .eq. 0) then
            cycle
        end if

        call as_mmhnme(idfimd, nomamd, edconn, edmail, tygeom, &
                       ednoda, nmatyp, codre2)
!
        if (codre2 .eq. 0) then
!
! 2.2.2. ==> SI LE CHOIX S'EST FAIT AVEC UNE VALEUR D'INSTANT, ON REPERE
!            LE NUMERO D'ORDRE ASSOCIE
!
            if (iinst .ne. 0) then
                if (niv .gt. 1) then
                    write (ifm, *) '.... INSTANT : ', inst
                end if
                call mdchin(nofimd, idfimd, nochmd, typent, tygeom, prefix, npas, codret)
                if (npas .ne. 0) then
                    call jeveuo(prefix//'.INST', 'L', vr=vinst)
                    call jeveuo(prefix//'.NUME', 'L', vi=nume)
                    logaux = ASTER_FALSE
                    do iaux2 = 1, npas
                        if (crit(1:4) .eq. 'RELA') then
                            if (abs(vinst(iaux2)-inst) .le. abs(prec*inst)) then
                                logaux = ASTER_TRUE
                            end if
                        else if (crit(1:4) .eq. 'ABSO') then
                            if (abs(vinst(iaux2)-inst) .le. abs(prec)) then
                                logaux = ASTER_TRUE
                            end if
                        end if
                        if (logaux) then
                            numpt = nume(1+2*iaux2-2)
                            numord = nume(1+2*iaux2-1)
                            goto 222
                        end if
                    end do
                    valk(1) = nofimd(1:24)
                    valk(2) = nochmd(1:24)
                    valr = inst
                    vali(1) = typent
                    vali(2) = typgeo(1)
                    call utmess('A', 'MED_97', nk=2, valk=valk, ni=2, vali=vali, sr=valr)
                    call utmess('A', 'MED_52')
                    goto 22
222                 continue
                    if (niv .gt. 1) then
                        valk(1) = nochmd(1:24)
                        vali(1) = typent
                        vali(2) = typgeo(1)
                        vali(3) = numord
                        vali(4) = numpt
                        valr = inst
                        call utmess('I', 'MED_86', sk=valk(1), ni=4, vali=vali, &
                                    sr=valr)
                    end if
                    call jedetr(prefix//'.INST')
                    call jedetr(prefix//'.NUME')
                end if
            end if
!
! 2.2.3. ==> RECHERCHE DES COMPOSANTES
!
            call codent(letype, 'G', k2bid)
            nmcmfl = '&&'//nompro//'.NOMCMP_FICHIE'//k2bid
!
            call mdexch(nofimd, idfimd, nochmd, numpt, numord, &
                        nbcmpv, ncmpvm, nbvato, typent, tygeom, &
                        existc, nbcmfi, nmcmfl, nbval, nbprof, &
                        codret)
            if (existc .ge. 3) then
                existt = ASTER_TRUE
                nbtylu = nbtylu+1
                nmcmfi(nbtylu) = nmcmfl
                if (typech(1:4) .ne. 'NOEU') then
                    lypent(nbtylu) = typent
                    lygeom(nbtylu) = tygeom
                    nlyval(nbtylu) = nbval
                    ltyp(nbtylu) = iaux
                    nbty(nbtylu) = nmatyp
                    lnbpro(nbtylu) = nbprof
                end if
            end if
!
!       ENDIF <<< IF ( CODRE2.EQ.0 )
        end if
!       INCREMENTE LE NUMERO INITIAL DES MAILLES DU TYPE SUIVANT
        if (nmatyp .gt. 0) then
            numma = numma+nmatyp
        end if
22      continue
    end do
!
! 2.3. ==> IL MANQUE DES CHOSES !
!
    if (.not. existt) then
        valk(1) = nofimd(1:64)
        valk(2) = nochmd(1:64)
        call utmess('A+', 'MED_98', nk=2, valk=valk)
        if (iinst .ne. 0) then
            valr = inst
            call utmess('A+', 'MED_68', sr=valr)
        else
            vali(1) = numord
            vali(2) = numpt
            call utmess('A+', 'MED_69', ni=2, vali=vali)
        end if
        if (existc .eq. 0) then
            call utmess('A', 'MED_32', sk=nochmd)
        else if (existc .eq. 1) then
            call utmess('A', 'MED_33')
        else if (existc .eq. 2) then
            if (iinst .ne. 0) then
                call utmess('A', 'MED_34')
            else
                call utmess('A', 'MED_35')
            end if
        else if (existc .eq. 4) then
            call utmess('A', 'MED_36')
        end if
        call utmess('F', 'MED_37')
    end if
!
!====
! 0. TRAITEMENT PARTICULIER POUR LES CHAMPS ELGA
!====
!
! 0.1 ==> VERIFICATION DES LOCALISATIONS DES PG
!         CREATION DU TABLEAU DEFINISSANT LE NBRE DE PG PAR MAIL
!         CREATION D'UN TABLEAU DE CORRESPONDANCE ENTRE LES PG MED/ASTER
!
    do i = 1, ntypel
        do j = 1, npgmax
            indpg(i, j) = 0
        end do
    end do
!
    if (typech(1:4) .eq. 'NOEU') then
        call cnscre(nomaas, nomgd, ncmprf, zk8(jnocmp), 'V', &
                    chames)
!
        call jeveuo(chames//'.CNSD', 'L', adsd)
        call jeveuo(chames//'.CNSV', 'E', adsv)
        call jeveuo(chames//'.CNSL', 'E', adsl)
    else
        if (typech(1:4) .eq. 'ELGA') then
            call lrmpga(nrofic, ligrel, nochmd, nbma, npgma, &
                        npgmm, nspmm, ntypel, npgmax, indpg, &
                        numpt, numord, option, param, nomaas)
            call cescre('V', chames, typech, nomaas, nomgd, &
                        ncmprf, zk8(jnocmp), npgma, nspmm, [-ncmprf])
        else if (typech(1:4) .eq. 'CART') then
            call cescre('V', chames, 'ELEM', nomaas, nomgd, &
                        ncmprf, zk8(jnocmp), [-1], [-1], [-ncmprf])
        else if (typech(1:4) .eq. 'ELNO' .or. typech(1:4) .eq. 'ELEM') then
            call cescre('V', chames, typech, nomaas, nomgd, &
                        ncmprf, zk8(jnocmp), [-1], [-1], [-ncmprf])
        else
            ASSERT(ASTER_FALSE)
        end if
!
        call jeveuo(chames//'.CESD', 'L', adsd)
        call jeveuo(chames//'.CESV', 'E', adsv)
        call jeveuo(chames//'.CESL', 'E', adsl)
    end if
!
!=====================================================================
! 3. TRAITEMENT DES CHAMPS AUX NOEUDS                             ====
!=====================================================================
!
    if (typech(1:4) .eq. 'NOEU') then
!
!====
! 3.1 LECTURE DES VALEURS
!====
!
        typent = typen
        tygeom = typnoe
        call lrcmle(idfimd, nochmd, nbcmfi, nbvato, numpt, &
                    numord, typent, tygeom, 1, ntvale, &
                    nomprf, codret)
!
!====
! 3.2   LECTURE DU PROFIL
!====
!
        if (nomprf .eq. ednopf) then
            lgproa = 0
        else
            call lrcmpr(idfimd, nomprf, ntproa, lgproa, codret)
        end if
!
!====
! 3.3   TRANFERT DES VALEURS
!====
!
        call lrcmva(ntvale, nbvato, ntproa, lgproa, ncmprf, &
                    zk8(jnocmp), nbcmfi, nmcmfi(1), nbcmpv, ncmpvm, &
                    numcmp, nochmd, adsl, adsv, codret)
!
    else
!
!=====================================================================
! 4. TRAITEMENT DES CHAMPS AUX ELEMENTS                           ====
!=====================================================================
!
!  ON BOUCLE SUR LES TYPES DE MAILLE LUES DANS LE CHAMP MED.
!  LES VALEURS NUMERIQUES SONT SAUVEES DANS LE TABLEAU D ADRESSE ADSV
!  CE TABLEAU A ETE DIMENSIONNE PAR CESCRE A :
!  NB DE TYPE DE MAIL * NB DE VALEURS PAR MAILLE * NB DE COMPOSANTES
!  * NB DE MAILLES DU TYPE
!  LE NB DE VALEURS PAR MAILLE :
!       - POUR UN ELNO : NB DE NOEUDS (NBNOMA DONNE PAR CONNECTIVITE)
!       - POUR UN ELEM : 1
!       - POUR UN ELGA : VARIABLE (INFO PRESENTE DANS LE TABLEAU NPGMA)
!
        ASSERT(zi(adsd) .eq. nbma)
        call jeveuo(nomaas(1:8)//'.TYPMAIL', 'L', vi=typmail)
!
        do letype = 1, nbtylu
!
            nbprof = lnbpro(letype)
            do iprof = 1, nbprof
                nbnoma = 1
                if (typech(1:4) .eq. 'ELNO') then
                    nbnoma = nnotyp(ltyp(letype))
                end if
                if (niv .gt. 1) then
                    write (ifm, *) '.... NBNOMA : ', nbnoma
                end if
!
!====
! 4.0   LECTURE DES VALEURS
!====
!
                call jedetr(ntvale)
                call lrcmle(idfimd, nochmd, nbcmfi, nlyval(letype), numpt, &
                            numord, lypent(letype), lygeom(letype), iprof, ntvale, &
                            nomprf, codret)
!
!====
! 4.1   LECTURE DU PROFIL
!====
!
                if (nomprf .eq. ednopf) then
                    lgproa = 0
                else
                    call jedetr(ntproa)
                    call lrcmpr(idfimd, nomprf, ntproa, lgproa, codret)
                    call jeveuo(ntproa, 'L', jntpro)
                    call jelira(ntproa, 'LONMAX', lgprof)
                end if
!
!====
! 4.2   VECTEUR CONTENANT LES NUMEROS DES MAILLES POUR CE TYPE
!====
!         ON BOUCLE (72) SUR LES MAILLES DU MAILLAGE ASTER
!         ET ON RELEVE LES MAILLES CORRESPONDANT AU TYPE LU
!
!         ON SOUHAITE VERIFIER QUE LE MODELE ASTER ET LE PROFIL
!         MED ONT BIEN LE MEME NOMBRE DE MAILLE DE CHAQUE TYPE
                call jeexin(ligrel//'.LGRF', iret)
                if (iret .ne. 0) then
                    call jeveuo(ligrel//'.LGRF', 'L', vk8=lgrf)
                    modele = lgrf(2)
                    call jeveuo(modele//'.MAILLE', 'L', jmaill)
                else
                    jmaill = 0
                end if
!
                call wkvect('&&'//nompro//'.NUM.'//nomtyp(ltyp(letype)), &
                            'V V I', nbty(letype), jnumty)
                k = 0
                if (lgproa .eq. 0) then
                    do ima = 1, nbma
                        if (typmail(ima) .eq. ltyp(letype)) then
                            if (jmaill .eq. 0) then
                                k = k+1
                                zi(jnumty+k-1) = ima
                            else if (zi(jmaill+ima-1) .ne. 0) then
                                k = k+1
                                zi(jnumty+k-1) = ima
                            end if
                        end if
                    end do
                    if (k .ne. nbty(letype)) then
                        call utmess('F', 'MED_58')
                    end if
                else
                    k = 0
                    cptyma = 1
                    do ima = 1, nbma
                        if (typmail(ima) .eq. ltyp(letype)) then
                            if (zi(jntpro+k) .eq. cptyma) then
                                if (jmaill .eq. 0) then
                                    k = k+1
                                    zi(jnumty+k-1) = ima
                                else if (zi(jmaill+ima-1) .ne. 0) then
                                    k = k+1
                                    zi(jnumty+k-1) = ima
                                end if
                            end if
                            cptyma = cptyma+1
                        end if
                    end do
                    if (k .ne. lgprof) then
                        call utmess('F', 'MED_58')
                    end if
                end if
!
!====
! 4.3   TRANFERT DES VALEURS
!====
!
                lrenum = ASTER_FALSE
                if (modnum(ltyp(letype)) .eq. 1) lrenum = ASTER_TRUE
                call lrcmve(ntvale, nbty(letype), nbnoma, ntproa, lgproa, &
                            ncmprf, zk8(jnocmp), ntypel, npgmax, indpg, &
                            nbcmfi, nmcmfi(letype), nbcmpv, ncmpvm, numcmp, &
                            jnumty, nochmd, nbma, npgma, npgmm, &
                            nspmm, typech, ltyp(letype), adsl, adsv, &
                            adsd, lrenum, nuanom, codret)
                call jedetr('&&'//nompro//'.NUM.'//nomtyp(ltyp(letype)))
!
            end do
        end do
    end if
!
!====
! 5. FIN
!====
!
! 5.1. ==> FERMETURE FICHIER
!
    call as_mficlo(idfimd, codret)
    if (codret .ne. 0) then
        saux08 = 'mficlo'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
!
! 5.2. ==> MENAGE
!
    call jedetr(numcmp)
    call jedetr(ntncmp)
    call jedetr(ntucmp)
    call jedetr(ntvale)
    call jedetr(ntproa)
!
    do letype = 0, nbtyp
        call codent(letype, 'G', k2bid)
        call jedetr('&&'//nompro//'.NOMCMP_FICHIE'//k2bid)
    end do
!
    if (niv .gt. 1) then
        write (ifm, 101) 'FIN DE '//nompro
    end if
    call jedema()
!
end subroutine
