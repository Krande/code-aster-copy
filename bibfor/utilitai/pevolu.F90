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

subroutine pevolu(resu, modele, carele, nbocc)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/chpchd.h"
#include "asterfort/chsut1.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnscno.h"
#include "asterfort/codent.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/exlim1.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nopar2.h"
#include "asterfort/pebpct.h"
#include "asterfort/reliem.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/rsutnu.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utflmd.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    integer(kind=8) :: nbocc
    character(len=8) :: modele, carele
    character(len=19) :: resu
!
!
!     OPERATEUR :  POST_ELEM
!     TRAITEMENT DU MOT CLE-FACTEUR : "VOLUMOGRAMME"
!     ------------------------------------------------------------------
!
    integer(kind=8) :: nr, nd, np, nc, ni, no, nli, nlo, iret, ibid, nbma, nbordr, jno
    integer(kind=8) :: nn, nbmaf, jma
    integer(kind=8) :: nbpar, nbpmax, iocc, inum, numo, jin, nbmato, iresma, ncmpm, ifm
    integer(kind=8) :: nbcmp, nbint, jbpct, ivalr, ii, i, ib, jvalr, jvali, jvalk, niv
    integer(kind=8) :: nucmp, ivali, bfix, ivol(2), tord(1)
    parameter(nbpmax=13)
    character(len=4) :: tych, ki
    character(len=8) :: mailla, crit, k8b, resuco, chamg, typpar(nbpmax), nomgd
    character(len=8) :: typmcl(1), tout, nomcmp, infoma, ncpini
    character(len=8) :: nopar, norme
    real(kind=8) :: prec, inst, borne(2), voltot, seuil
    complex(kind=8) :: c16b
    character(len=19) :: knum, kins, lisins, cham, cham2, chamtm, celmod, ligrel
    character(len=19) :: tmpcha, cham3, cespoi
    character(len=16) :: nompar(nbpmax), mocles(1), optio2, nomcha, valk, valr
    character(len=16) :: vali
    character(len=24) :: mesmai, mesmaf, mesmae, borpct, valk2(5), grouma
    aster_logical :: exiord, toneut, lseuil
    character(len=8), pointer :: cmp1(:) => null()
    character(len=8), pointer :: cmp2(:) => null()
    character(len=8), pointer :: cnsc(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
!
! --- 1- RECUPERATION DU MAILLAGE ET DU NOMBRE DE MAILLES
!     ===================================================
    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mailla)
    call dismoi('NB_MA_MAILLA', mailla, 'MAILLAGE', repi=nbmato)
!
!
! --- 2- RECUPERATION DU RESULTAT ET DES NUMEROS D'ORDRE
!     ==================================================
    call getvid(' ', 'RESULTAT', scal=resuco, nbret=nr)
    call getvid(' ', 'CHAM_GD', scal=chamg, nbret=nd)
!
    call getvr8(' ', 'PRECISION', scal=prec, nbret=np)
    call getvtx(' ', 'CRITERE', scal=crit, nbret=nc)
    call getvr8(' ', 'INST', nbval=0, nbret=ni)
    call getvis(' ', 'NUME_ORDRE', nbval=0, nbret=no)
    call getvid(' ', 'LIST_INST', nbval=0, nbret=nli)
    call getvid(' ', 'LIST_ORDRE', nbval=0, nbret=nlo)
!
    valr = '&&PEVOLU.VALR'
    vali = '&&PEVOLU.VALI'
    valk = '&&PEVOLU.VALK'
    knum = '&&PEVOLU.NUME_ORDRE'
    kins = '&&PEVOLU.INST'
    mesmai = '&&PEVOLU.MES_MAILLES'
    ligrel = '&&PEVOLU.LIGREL'
    mesmaf = '&&PEVOLU.MAILLES_FILTRE'
    borpct = '&&PEVOLU_BORNES_PCT'
    cespoi = '&&PEVOLU_POIDS'
!
    exiord = .false.
!
    if (nd .ne. 0) then
!
        nbordr = 1
        call wkvect(knum, 'V V I', nbordr, jno)
        zi(jno) = 1
        exiord = .true.
!
    else
!
!       -- NUME_ORDRE --
        if (no .ne. 0) then
            exiord = .true.
            nbordr = -no
            call wkvect(knum, 'V V I', nbordr, jno)
            call getvis(' ', 'NUME_ORDRE', nbval=nbordr, vect=zi(jno), nbret=iret)
        end if
!
!       -- LIST_ORDRE --
        if (nlo .ne. 0) then
            exiord = .true.
            call getvid(' ', 'LIST_ORDRE', scal=lisins, nbret=iret)
            call jeveuo(lisins//'.VALE', 'L', jno)
            call jelira(lisins//'.VALE', 'LONMAX', nbordr)
        end if
!
!       -- INST --
        if (ni .ne. 0) then
            nbordr = -ni
            call wkvect(kins, 'V V R', nbordr, jin)
            call getvr8(' ', 'INST', nbval=nbordr, vect=zr(jin), nbret=iret)
        end if
!
!       -- LIST_INST --
        if (nli .ne. 0) then
            call getvid(' ', 'LIST_INST', scal=lisins, nbret=iret)
            call jeveuo(lisins//'.VALE', 'L', jin)
            call jelira(lisins//'.VALE', 'LONMAX', nbordr)
        end if
!
!       -- TOUT_ORDRE --
        nn = nli+ni+no+nlo
        if (nn .eq. 0) then
            exiord = .true.
            call rsutnu(resuco, ' ', 0, knum, nbordr, &
                        prec, crit, iret)
            call jeveuo(knum, 'L', jno)
        end if
!
    end if
!
!
! --- 3- CREATION DE LA TABLE
!     =======================
    call tbcrsd(resu, 'G')
    if (nr .ne. 0) then
        nbpar = 11
        nompar(1) = 'RESULTAT'
        nompar(2) = 'NOM_CHAM'
        nompar(3) = 'NUME_ORDRE'
        nompar(4) = 'INST'
        nompar(5) = 'NOM_CMP'
        nompar(6) = 'GROUP_MA'
        nompar(7) = 'RESTRICTION'
        nompar(8) = 'INTERVALLE'
        nompar(9) = 'BORNE_INF'
        nompar(10) = 'BORNE_SUP'
        nompar(11) = 'DISTRIBUTION'
        typpar(1) = 'K8'
        typpar(2) = 'K16'
        typpar(3) = 'I'
        typpar(4) = 'R'
        typpar(5) = 'K8'
        typpar(6) = 'K24'
        typpar(7) = 'K8'
        typpar(8) = 'I'
        typpar(9) = 'R'
        typpar(10) = 'R'
        typpar(11) = 'R'
    else
        nbpar = 8
        nompar(1) = 'CHAM_GD'
        nompar(2) = 'NOM_CMP'
        nompar(3) = 'GROUP_MA'
        nompar(4) = 'RESTRICTION'
        nompar(5) = 'INTERVALLE'
        nompar(6) = 'BORNE_INF'
        nompar(7) = 'BORNE_SUP'
        nompar(8) = 'DISTRIBUTION'
        typpar(1) = 'K8'
        typpar(2) = 'K8'
        typpar(3) = 'K24'
        typpar(4) = 'K8'
        typpar(5) = 'I'
        typpar(6) = 'R'
        typpar(7) = 'R'
        typpar(8) = 'R'
    end if
    call tbajpa(resu, nbpar, nompar, typpar)
!
! --- 4- REMPLISSAGE DE LA TABLE
!     ==========================
!
!     --- ON PARCOURT LES OCCURENCES DU MOT CLE 'VOLUMOGRAMME':
!     ---------------------------------------------------------
    if (nr .eq. 0) then
        tmpcha = 'TMP_CHAMP_GD'
        call copisd('CHAMP', 'V', chamg, tmpcha)
    end if
!
    do iocc = 1, nbocc

!       -- 4.1 recuperation des mailles, calcul de ligrel :
!       -----------------------------------------------------

        call getvtx('VOLUMOGRAMME', 'TOUT', iocc=iocc, scal=tout, nbret=iret)
        if (iret .ne. 0) then
            mocles(1) = 'TOUT'
            typmcl(1) = 'TOUT'
            grouma = '-'
        else
            mocles(1) = 'GROUP_MA'
            typmcl(1) = 'GROUP_MA'
            call getvtx('VOLUMOGRAMME', 'GROUP_MA', iocc=iocc, scal=grouma, nbret=iret)
        end if
!
!       -- MAILLES FOURNIES PAR L'UTILISATEUR -
        call reliem(modele, mailla, 'NU_MAILLE', 'VOLUMOGRAMME', iocc, &
                    1, mocles, typmcl, mesmai, nbma)
!
        mesmae = mesmai
!       -- MAILLES EVENTUELLEMENT FILTREES EN FONCTION DE LA DIMENSION
!       GEOMETRIQUE (2D OU 3D)
        call getvtx('VOLUMOGRAMME', 'TYPE_MAILLE', iocc=iocc, scal=infoma, nbret=iret)
        if (iret .ne. 0) then
            if (infoma .eq. '2D') then
                iresma = 2
            else if (infoma .eq. '3D') then
                iresma = 3
            else
                ASSERT(.false.)
            end if
            call utflmd(mailla, mesmai, nbma, iresma, ' ', &
                        nbmaf, mesmaf)
            if (nbmaf .gt. 0) then
                call jedetr(mesmai)
                nbma = nbmaf
                mesmae = mesmaf
            else
                call utmess('F', 'PREPOST2_6')
            end if
        else
            infoma = '-'
        end if

!       -- calcul de ligrel :
!       ---------------------
        call jeveuo(mesmae, 'L', jma)
        call jelira(mesmae, 'LONMAX', nbma)
        call exlim1(zi(jma), nbma, modele, 'V', ligrel)

        if (nr .ne. 0) then
            call getvtx('VOLUMOGRAMME', 'NOM_CHAM', iocc=iocc, scal=nomcha, nbret=iret)
            if (iret .eq. 0) then
                call utmess('F', 'POSTELEM_4')
            end if
        else
            nomcha = chamg
            cham2 = tmpcha
        end if

!       -- BOUCLE SUR LES NUMEROS D'ORDRE:
!       ----------------------------------
!
        do inum = 1, nbordr
            toneut = .false.

!
!           -- 4.2 RECUPERATION DU CHAMP --
!
            if (nr .ne. 0) then
!               --  RESULTAT --
                if (exiord) then
!                  - ORDRE -
                    numo = zi(jno+inum-1)
                    call rsadpa(resuco, 'L', 1, 'INST', numo, &
                                0, sjv=jin, styp=k8b)
                    inst = zr(jin)
                else
!                   - INST -
                    inst = zr(jin+inum-1)
                    call rsorac(resuco, 'INST', 0, zr(jin+inum-1), k8b, &
                                c16b, prec, crit, tord, nbordr, &
                                iret)
                    numo = tord(1)
                end if
!
                call rsexch('F', resuco, nomcha, numo, cham2, &
                            iret)
!
            else
!               -- CHAM_GD --
                numo = nbordr
            end if
!
            call dismoi('TYPE_CHAMP', cham2, 'CHAMP', repk=tych, arret='C', &
                        ier=iret)
            call dismoi('NOM_GD', cham2, 'CHAMP', repk=nomgd, arret='C', &
                        ier=iret)
!
            if (nomgd(6:6) .eq. 'C') goto 10
!
            if (tych(1:2) .ne. 'EL') then
!
!          --- 1. TRANSFORMATION DU CHAMP EN CHAMP NEUTRE:
!              - CHANGEMENT DE LA GRANDEUR EN NEUT_R
!              - CHAMGEMENT DES COMPOSANTES EN X1,X2,X3,...
                toneut = .true.
                chamtm = '&&PEVOLU.CHS1'
                call cnocns(cham2, 'V', chamtm)
                call jeveuo(chamtm//'.CNSC', 'L', vk8=cnsc)
                call jelira(chamtm//'.CNSC', 'LONMAX', ncmpm)
                AS_ALLOCATE(vk8=cmp1, size=ncmpm)
                AS_ALLOCATE(vk8=cmp2, size=ncmpm)
                do i = 1, ncmpm
                    call codent(i, 'G', ki)
                    cmp2(i) = 'X'//ki(1:len(ki))
                    cmp1(i) = cnsc(i)
                end do
                call chsut1(chamtm, 'NEUT_R', ncmpm, cmp1, cmp2, &
                            'V', chamtm)
!
                cham3 = '&&PEEINT.CHAM_3'
                call cnscno(chamtm, ' ', 'NON', 'V', cham3, &
                            'F', ibid)
                call detrsd('CHAM_NO_S', chamtm)
!
!           --- 2. CHANGEMENT DE DISCRETISATION : NOEU -> ELGA
                optio2 = 'TOU_INI_ELGA'
                call dismoi('NOM_GD', cham3, 'CHAMP', repk=nomgd, arret='C', &
                            ier=iret)
                call nopar2(optio2, nomgd, 'OUT', nopar)
                celmod = '&&PEVOLU.CELMOD'
                call alchml(ligrel, optio2, nopar, 'V', celmod, &
                            ib, ' ')
                if (ib .ne. 0) then
                    valk2(1) = ligrel
                    valk2(2) = nopar
                    valk2(3) = optio2
                    call utmess('F', 'UTILITAI3_23', nk=3, valk=valk2)
                end if
                cham = '&&CHPCHD.CHAM'
                call chpchd(cham3, 'ELGA', celmod, 'OUI', 'V', cham, modele)
                call detrsd('CHAMP', celmod)
                call detrsd('CHAMP', cham3)
!
            else
                cham = cham2
            end if
!
            call dismoi('TYPE_CHAMP', cham, 'CHAMP', repk=tych, arret='C', &
                        ier=iret)

!
!      -- 4.3 RECUPERATION DE LA COMPOSANTE --
!
            call getvtx('VOLUMOGRAMME', 'NOM_CMP', iocc=iocc, scal=nomcmp, nbret=nbcmp)
            ncpini = nomcmp
            if (toneut) then
                nucmp = indik8(cmp1, nomcmp, 1, ncmpm)
                nomcmp = cmp2(nucmp)
                AS_DEALLOCATE(vk8=cmp1)
                AS_DEALLOCATE(vk8=cmp2)
            end if
!
!      -- 4.4 CALCUL DES INTERVALLES ET DE LA DISTRIBUTION --
!
!        - ON RECUPERE LE NOMBRE D'INTERVALLES OU LE SEUIL ET LA NORME
!        - POUR CALCUL RELATIF OU ABSOLU
            call getvis('VOLUMOGRAMME', 'NB_INTERV', iocc=iocc, scal=nbint, nbret=iret)
            seuil = 0.d0
            lseuil = .false.
            call getvr8('VOLUMOGRAMME', 'SEUIL', iocc=iocc, scal=seuil, nbret=iret)
            if (iret .ne. 0) then
                nbint = 2
                lseuil = .true.
            end if
            call getvtx('VOLUMOGRAMME', 'NORME', iocc=iocc, scal=norme, nbret=iret)
            call getvr8('VOLUMOGRAMME', 'BORNES', iocc=iocc, nbval=2, vect=borne, &
                        nbret=iret)
            if (iret .ne. 0) then
                bfix = 1
                ASSERT(borne(1) .lt. borne(2))
            else
                borne(1) = 0.d0
                borne(2) = 0.d0
                bfix = 0
            end if
!        - ON CREE UN TABLEAU POUR RECUEILLIR POUR CHAQUE INTERVALLE
!         LES 3 VALEURS SUIVANTES:
!             . LA BORNE INF,
!             . LA BORNE SUP,
!             . LE VALEUR DE REPARTITION DE LA COMPOSANTE
            call wkvect(borpct, 'V V R', 3*nbint, jbpct)
            call pebpct(ligrel, nbma, mesmae, cham, nomcmp, &
                        3*nbint, bfix, borne, norme, seuil, &
                        lseuil, zr(jbpct), voltot, carele, cespoi)
!
!      -- 4.5 ON REMPLIT LA TABLE --
!
            if (nompar(1) .eq. 'RESULTAT') then
                call wkvect(valk, 'V V K24', 5, jvalk)
                zk24(jvalk) = resuco
                zk24(jvalk+1) = nomcha
                zk24(jvalk+2) = ncpini
                zk24(jvalk+3) = grouma
                zk24(jvalk+4) = infoma
                call wkvect(valr, 'V V R', 4, jvalr)
                zr(jvalr) = inst
                ivalr = 1
                call wkvect(vali, 'V V I', 2, jvali)
                zi(jvali) = numo
                ivali = 1
            else
                call wkvect(valk, 'V V K24', 4, jvalk)
                zk24(jvalk) = nomcha
                zk24(jvalk+1) = ncpini
                zk24(jvalk+2) = grouma
                zk24(jvalk+3) = infoma
                call wkvect(valr, 'V V R', 3, jvalr)
                ivalr = 0
                call wkvect(vali, 'V V I', 1, jvali)
                ivali = 0
            end if
!
!        - POUR CHAQUE INTERVALLE, ON AJOUTE UNE LIGNE A LA TABLE :
            do ii = 1, nbint
                zr(jvalr+ivalr) = zr(jbpct+3*(ii-1))
                zr(jvalr+ivalr+1) = zr(jbpct+3*(ii-1)+1)
                zr(jvalr+ivalr+2) = zr(jbpct+3*(ii-1)+2)
                zi(jvali+ivali) = ii
                call tbajli(resu, nbpar, nompar, zi(jvali), zr(jvalr), &
                            [c16b], zk24(jvalk), 0)
            end do
!
!     IMPRESSION DU VOLUME TOTAL CONCERNE PAR LE CALCUL :
!      - SOIT LE VOLUME TOTAL DEFINI PAR L ENTITE TOPOLOGIQUE
!      - SOIT POUR CELUI DONT LES VALEURS SONT COMPRISES DANS LES
!        BORNES SI CELLES-CI SONT RENSEIGNEES PAR L UTILISATEUR
            ivol(1) = iocc
            ivol(2) = inum
            call utmess('I', 'UTILITAI7_14', ni=2, vali=ivol, sr=voltot)
!
!      -- 4.5 NETTOYAGE POUR L'OCCURRENCE SUIVANTE --
            call jedetr(valr)
            call jedetr(vali)
            call jedetr(valk)
            call jedetr(borpct)
!
!     --- FIN DE LA BOUCLE SUR LES NUMEROS D'ORDRE:
!     ---------------------------------------------
        end do

        call jedetr(mesmaf)
        call jedetr(mesmai)
        call detrsd('LIGREL', ligrel)
        call detrsd('CHAM_ELEM_S', cespoi)

!     --- FIN DE LA BOUCLE SUR LES OCCURRENCES DU MOT-CLE VOLUMOGRAMME
!     ----------------------------------------------------------------
10      continue
    end do
!
    if (nr .eq. 0) then
        call detrsd('CHAMP', tmpcha)
    end if
!
!
    call jedema()
!
end subroutine
