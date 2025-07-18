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
subroutine op0151()
    implicit none
!       C A L C _ F A T I G U E
!       ----------------------------------------------------------------
!       CREATION D UN CHAM_ELEM D ISODOMMAGE A LA FATIGUE
!       D UN MATERIAU SOUMIS A UN CYCLAGE EN CONTRAINTES
!       A PARTIR D'UN CHAM_ELEM DE GRANDEUR 1D EQUIVALENTE
!       ----------------------------------------------------------------
!       ATTENTION : LE CHAM_ELEM EN SORTIE EST CONSTRUIT ARTIFICIELEMENT
!                   SANS PASSER PAR LA ROUTINE CALCUL
!                   A PARTIR DU CHAM_ELEM DANS LA SD RESULTAT
!                   DE TYPE EVOL_NOLI , EVOL_ELAS , DYNA_TRANS, CREE
!                   PAR CALC_CHAMP (OPTIONS SIEQ(EPEQ)_ELNO(ELGA))
!                   LA COHERENCE DU CHAMP EST DONC LIEE A CE DERNIER
!                   (MEMES ELEMENTS/GREL POUR LE CALCUL DE L OPTION ...)
!
!       IMPLANTE  : ACTUELLEMENT  :
!                  - CALCUL DU DOMMAGE A PARTIR DE   = /CONTRAINTE
!                                                      /DEFORMATION
!                  - POINTS DE CALCUL DU DOMMAGE     = /NOEUDS
!                                                      /POINTS DE GAUSS
!                  - COMPOSANTES GRANDEUR EQUIVALENTE= VMIS_SG....
!                  - METHODE  D'EXTRACTION DES PICS  = RAINFLOW
!                  - METHODE  DE COMPTAGE DES CYCLES = RAINFLOW
!                  - METHODE  DE CALCUL   DU DOMMAGE = /WOHLER
!                                                      /MANSON_COFFIN
!                                                      /TAHERI_MANSON
!                                                      /TAHERI_MIXTE
!                  - OPTIONS OUT                     = /DOMA_ELNO_SIGM
!                                                      /DOMA_ELGA_SIGM
!                                                      /DOMA_ELGA_EPSI
!                                                      /DOMA_ELNO_EPSI
!                                                      /DOMA_ELGA_EPME
!                                                      /DOMA_ELNO_EPME
!       ----------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/alchml.h"
#include "asterfort/anacri.h"
#include "asterfort/assert.h"
#include "asterfort/dmgmod.h"
#include "asterfort/fgvdmg.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeimpo.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/paqmai.h"
#include "asterfort/paqnoe.h"
#include "asterfort/rccome.h"
#include "asterfort/rcpare.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
!       ---------------------------------------------------------------
    integer(kind=8) :: icodre, icodwo, icodba, icodhs, icodma
    character(len=8) :: nomu, nomres, nommai, nommat
    character(len=8) :: nomfon, nomnap, nommod, inscri
    character(len=16) :: concep, cmd, pheno, typcal, nomcri, nommet, cara
    character(len=16) :: proaxe, nomsym, typcha, nomopt, nomgde, nomfor, forcri
    character(len=16) :: forvie, grdvie
    character(len=16) :: mexpic, mcompt, mdomag, typeq, typoi, typdg, option
    character(len=19) :: nomsd, chelem, chelrs, ligrel, nomsd2
    character(len=24) :: valk(6)
    aster_logical :: fordef, crsigm, crepst, crepse, crepsp
    real(kind=8) :: prec, instic
!
    integer(kind=8) :: nval, impr, ifm, jordr, jcoef, jcelk
    integer(kind=8) :: nbpt, nbord, nbcmp, numcmp(6), ntcmp, ibid
    integer(kind=8) :: ivdmg, numsym, nbpt2, nbord2, iret, ivch
    integer(kind=8) :: vali(2), paract(35)
!
! ------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i
    real(kind=8), pointer :: celv(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    call infmaj()
!
    nomfon = ' '
    nomnap = ' '
!
! ----- DONNEES UTILISATEUR
!
    call getres(nomu, concep, cmd)
!
! ---   TYPE DE CALCUL
    call getvtx(' ', 'TYPE_CALCUL', scal=typcal, nbret=nval)
!
! ---------------------------------------------------------------------
! ---- FATIGUE MULTIAXIALE
! ---------------------------------------------------------------------
!
    if (typcal(1:13) .eq. 'FATIGUE_MULTI') then
!
! ---   TYPE DU CHARGEMENT APPLIQUE (PERIODIQUE OU NON_PERIODIQUE)
!
        call getvtx(' ', 'TYPE_CHARGE', scal=typcha, nbret=nval)
!
! ---   NOM DE L'OPTION (CALCUL AUX POINTS DE GAUSS OU AUX NOEUDS
!                        CHAM_NO)
        call getvtx(' ', 'OPTION', scal=nomopt, nbret=nval)
!
! ---   STRUCTURE RESULTAT CONTENANT LES CHAM_ELEMS DE SIGMA
!       OU LES CHAM_NOS DE SIGMA
        call getvid(' ', 'RESULTAT', scal=nomres, nbret=nval)
        nomsd = nomres
!
! ---   NOM DU CRITERE
        call getvtx(' ', 'CRITERE', scal=nomcri, nbret=nval)
!
        call getvid(' ', 'FORMULE_GRDEQ', scal=nomfor, nbret=nval)
        if (nval .eq. 0) then
            nomfor = '        '
        end if
!
        call getvtx(' ', 'COURBE_GRD_VIE', scal=grdvie, nbret=nval)
        if (nval .eq. 0) then
            grdvie = '        '
        end if
!
        call getvid(' ', 'FORMULE_VIE', scal=forvie, nbret=nval)
        if (nval .eq. 0) then
            forvie = '        '
        end if
!
        call getvid(' ', 'FORMULE_CRITIQUE', scal=forcri, nbret=nval)
        if (nval .eq. 0) then
            forcri = '        '
        end if
!
!
! ---   NOM DE LA METHODE PERMETTANT DE DETERMINER LE CERCLE CIRCONSCRIT
        call getvtx(' ', 'METHODE', scal=nommet, nbret=nval)
        if (nval .eq. 0) then
            nommet = '        '
        end if
!
! ---   PROJECTION SUR UN AXE OU SUR DEUX AXES
!       (CHARGEMENT NON_PERIODIQUE UNIQUEMENT)
        call getvtx(' ', 'PROJECTION', scal=proaxe, nbret=nval)
        if (nval .eq. 0) then
            proaxe = '        '
        end if
!
! ---   NOM DU MAILLAGE
        call getvid(' ', 'MAILLAGE', scal=nommai, nbret=nval)
        if (nval .eq. 0) then
            nommai = '        '
        end if
!
! ----- INSTANT INITIAL DU PARTIE DE CHARGEMENT CYCLIQUE
!
        prec = 1.d-6
        inscri = '        '
        instic = 0.d0
!
        call getvr8(' ', 'INST_INIT_CYCL', scal=instic, nbret=nval)
        call getvtx(' ', 'INST_CRIT', scal=inscri, nbret=nval)
        call getvr8(' ', 'PRECISION', scal=prec, nbret=nval)
!
!
!---    ANALYSER LE CRITERE
! C  INITIALISER
        crsigm = .false.
        crepst = .false.
        crepse = .false.
        crepsp = .false.
        call anacri(nomcri, nomfor, typcha, 'OUI', paract, &
                    fordef, crsigm, crepst, crepse, crepsp)
!
!   FORDEF EST UNE BOOLEAN QUI INDIQUE S'IL EXISTE LE PARAMETRE
!   DE DEFORMATION DAS LA FORMULE (COMME DANS FATEMISOCIE)
!
        if (nomopt .eq. 'DOMA_ELGA') then
!
! ---   CONSTRUCTION DES PAQUETS DE MAILLES
            call paqmai(nomres, nomu, nommai, nommet, nomcri, &
                        nomfor, grdvie, forvie, forcri, fordef, &
                        typcha, proaxe, instic, inscri, prec)
!
        else if (nomopt .eq. 'DOMA_NOEUD') then
!
! ---   CONSTRUCTION DES PAQUETS DE NOEUDS
            call paqnoe(nomres, nomu, nommai, nommet, nomcri, &
                        nomfor, grdvie, forvie, forcri, fordef, &
                        typcha, proaxe, instic, inscri, prec)
        end if
!
!
        goto 7777
    end if
!
! ---------------------------------------------------------------------
! ---- CAS GENERAL (CUMUL DE DOMMAGE OU FATIGUE MODALE)
! ---------------------------------------------------------------------
!
! ---   NOM DE LA GRANDEUR EQUIVALENTE
    call getvtx('HISTOIRE', 'EQUI_GD', iocc=1, scal=nomgde, nbret=nval)
!
! ---   IMPRESSIONS
    call getvis(' ', 'INFO', scal=impr, nbret=nval)
!
! ---   CHAMP : NOM DE L'OPTION RESULTANTE
!   'DOMA_ELNO_SIGM'/'DOMA_ELGA_SIGM'/'DOMA_ELNO_EPSI'/'DOMA_ELGA_EPSI
!   'DOMA_ELNO_EPME'/'DOMA_ELGA_EPME'
!
    call getvtx(' ', 'OPTION', scal=nomopt, nbret=nval)
!
! ---   NOM DE LA METHODE DE CALCUL DU DOMMAGE
    call getvtx(' ', 'DOMMAGE', scal=mdomag, nbret=nval)
    call getvid(' ', 'MATER', scal=nommat, nbret=nval)
!
! ---   NOMBRE DE NUMEROS D ORDRE DE POINTS (NOEUDS/PG) DE CMPS ...
!       ---------------------------------------------------------------
!       NOM SYMBOLIQUE DE L OPTION IN UTILISEE (GRANDEUR EQUIVALENTE)
!       (SIEQ_ELNO  SIEQ_ELGA  EPEQ_ELGA EPEQ_ELNO
!        EPMQ_ELNO  EPMQ_ELGA  )
!       ET NOMBRE TOTAL DE COMPOSANTES DE CETTE OPTION
!
    if (nomopt(11:14) .eq. 'SIGM') then
        nomsym = 'SIEQ_'//nomopt(6:9)
    else if (nomopt(11:14) .eq. 'EPSI') then
        nomsym = 'EPEQ_'//nomopt(6:9)
    else if (nomopt(11:14) .eq. 'EPME') then
        nomsym = 'EPMQ_'//nomopt(6:9)
    else
        ASSERT(.false.)
    end if
!
!
    if (nomopt(6:14) .eq. 'ELGA_SIGM') then
        ntcmp = 17
    else if (nomopt(6:14) .eq. 'ELNO_SIGM') then
        ntcmp = 17
    else if (nomopt(6:14) .eq. 'ELGA_EPSI') then
        ntcmp = 14
    else if (nomopt(6:14) .eq. 'ELNO_EPSI') then
        ntcmp = 14
    else if (nomopt(6:14) .eq. 'ELGA_EPME') then
        ntcmp = 14
    else if (nomopt(6:14) .eq. 'ELNO_EPME') then
        ntcmp = 14
    else
        ASSERT(.false.)
    end if
!
!       TYPE DE GRANDEUR EQUIVALENTE UTILISEE LE POUR CALCUL DU DOMMAGE
!       ET NOMBRE DE COMPOSANTES DE CETTE GRANDEUR
!       VMIS_SIG   = CMP NUMERO 6       DE  L OPTION EQUI_...._SIGM
!       INVA_2_SG  = CMP NUMERO 5       DE  L OPTION EQUI_...._EPSI
!       INVA_2_SG  = CMP NUMERO 5       DE  L OPTION EQUI_...._EPME
!
    if (nomgde(1:7) .eq. 'VMIS_SG') then
        numcmp(1) = 6
        nbcmp = 1
    else if (nomgde(1:9) .eq. 'INVA_2_SG') then
        numcmp(1) = 5
        nbcmp = 1
    end if
!
! ---------------------------------------------------------------------
! ---- CUMUL DE DOMMAGE
! ---------------------------------------------------------------------
!
    if (typcal(1:13) .eq. 'CUMUL_DOMMAGE') then
!
! ---   STRUCTURE RESULTAT CONTENANT LES CHAM_ELEMS DE SIGMA EQUIVALENT
        call getvid('HISTOIRE', 'RESULTAT', iocc=1, scal=nomres, nbret=nval)
        nomsd = nomres
!
        if (mdomag .eq. 'WOHLER') then
            if (nomopt(11:14) .ne. 'SIGM') then
                call utmess('F', 'FATIGUE1_29', sk=nomopt)
            end if
            pheno = 'FATIGUE'
            call rccome(nommat, pheno, icodre)
            if (icodre .eq. 1) then
                call utmess('F', 'FATIGUE1_24')
            end if
            cara = 'WOHLER'
            call rcpare(nommat, pheno, cara, icodwo)
            cara = 'A_BASQUIN'
            call rcpare(nommat, pheno, cara, icodba)
            cara = 'A0'
            call rcpare(nommat, pheno, cara, icodhs)
            if (icodwo .ne. 0 .and. icodba .ne. 0 .and. icodhs .ne. 0) then
                call utmess('F', 'FATIGUE1_30')
            end if
!
        else if (mdomag .eq. 'MANSON_COFFIN') then
            if (nomopt(11:14) .ne. 'EPSI' .and. nomopt(11:14) .ne. 'EPME') then
                call utmess('F', 'FATIGUE1_31', sk=nomopt)
            end if
            pheno = 'FATIGUE'
            call rccome(nommat, pheno, icodre)
            if (icodre .eq. 1) then
                call utmess('F', 'FATIGUE1_24')
            end if
            cara = 'MANSON_COFFIN'
            call rcpare(nommat, pheno, cara, icodma)
            if (icodma .ne. 0) then
                call utmess('F', 'FATIGUE1_32')
            end if
!
        else if (mdomag .eq. 'TAHERI_MANSON') then
            if (nomopt(11:14) .ne. 'EPSI' .and. nomopt(11:14) .ne. 'EPME') then
                call utmess('F', 'FATIGUE1_25', sk=nomopt)
            end if
            pheno = 'FATIGUE'
            call rccome(nommat, pheno, icodre)
            if (icodre .eq. 1) then
                call utmess('F', 'FATIGUE1_24')
            end if
            cara = 'MANSON_COFFIN'
            call rcpare(nommat, pheno, cara, icodma)
            if (icodma .ne. 0) then
                call utmess('F', 'FATIGUE1_32')
            end if
            call getvid(' ', 'TAHERI_NAPPE', scal=nomnap, nbret=nval)
            if (nval .eq. 0) then
                call utmess('F', 'FATIGUE1_26')
            end if
            call getvid(' ', 'TAHERI_FONC', scal=nomfon, nbret=nval)
            if (nval .eq. 0) then
                call utmess('F', 'FATIGUE1_27')
            end if
!
        else if (mdomag .eq. 'TAHERI_MIXTE') then
            if (nomopt(11:14) .ne. 'EPSI' .and. nomopt(11:14) .ne. 'EPME') then
                call utmess('F', 'FATIGUE1_28', sk=nomopt)
            end if
            pheno = 'FATIGUE'
            call rccome(nommat, pheno, icodre)
            if (icodre .eq. 1) then
                call utmess('F', 'FATIGUE1_24')
            end if
            cara = 'MANSON_COFFIN'
            call rcpare(nommat, pheno, cara, icodma)
            if (icodma .ne. 0) then
                call utmess('F', 'FATIGUE1_32')
            end if
            cara = 'WOHLER'
            call rcpare(nommat, pheno, cara, icodwo)
            cara = 'A_BASQUIN'
            call rcpare(nommat, pheno, cara, icodba)
            cara = 'A0'
            call rcpare(nommat, pheno, cara, icodhs)
            if (icodwo .ne. 0 .and. icodba .ne. 0 .and. icodhs .ne. 0) then
                call utmess('F', 'FATIGUE1_30')
            end if
            call getvid(' ', 'TAHERI_NAPPE', scal=nomnap, nbret=nval)
            if (nval .eq. 0) then
                call utmess('F', 'FATIGUE1_26')
            end if
!
        end if
!
! --- VERIFICATION DU NOMBRE DE PAS DE TEMPS
        call jelira(nomsd//'.ORDR', 'LONUTI', nbord)
        if (nbord .lt. 2) then
            call utmess('F', 'FATIGUE1_76', si=nbord)
        end if
!
        call jenonu(jexnom(nomsd//'.DESC', nomsym), numsym)
        if (numsym .eq. 0) then
            valk(1) = nomsym
            valk(2) = nomsd
            call utmess('F', 'PREPOST4_5', nk=2, valk=valk)
        end if
        call jeveuo(jexnum(nomsd//'.TACH', numsym), 'L', ivch)
        chelrs = zk24(ivch) (1:19)
        if (chelrs .eq. ' ') then
            valk(1) = chelrs
            valk(2) = nomsym
            valk(3) = nomsd
            call utmess('F', 'PREPOST4_6', nk=3, valk=valk)
        end if
        call jeveuo(chelrs//'.CELK', 'L', jcelk)
        ligrel = zk24(jcelk-1+1) (1:19)
        call jelira(chelrs//'.CELV', 'LONMAX', nval)
!
!  -      IL Y A NTCMP COMPOSANTES DANS L OPTION XXXX_EQUI_YYYY
        nbpt = nval/ntcmp
!
        if (impr .ge. 2) then
            vali(1) = nbord
            vali(2) = nbpt
            call utmess('I', 'PREPOST6_27', ni=2, vali=vali)
        end if
!
! ----- CALCUL DU VECTEUR DOMMAGE EN CHAQUE NOEUD/PG
!       ----------------------------------------------------------------
!
        call wkvect('&&OP0151.DOMMAGE', 'V V R', nbpt, ivdmg)
!
        mexpic = 'RAINFLOW'
        mcompt = 'RAINFLOW'
        if (mdomag(1:6) .eq. 'TAHERI') mcompt = 'TAHERI'
!
        if (impr .ge. 2) then
            typeq = nomgde
            if (nomopt(11:14) .eq. 'SIGM') typdg = 'CONTRAINTE'
            if (nomopt(11:14) .eq. 'EPSI') typdg = 'DEFORMATION'
            if (nomopt(11:14) .eq. 'EPME') typdg = 'DEFORMATION'
            if (nomopt(6:9) .eq. 'ELNO') typoi = 'NOEUDS'
            if (nomopt(6:9) .eq. 'ELGA') typoi = 'POINTS DE GAUSS'
            valk(1) = typdg
            valk(2) = typoi
            valk(3) = typeq
            valk(4) = mexpic
            valk(5) = mcompt
            valk(6) = mdomag
            call utmess('I', 'PREPOST6_28', nk=6, valk=valk)
        end if
!
        call fgvdmg(nomsym, nomsd, nommat, nomnap, nomfon, &
                    mexpic, mcompt, mdomag, nbord, nbpt, &
                    ntcmp, nbcmp, numcmp, impr, zr(ivdmg))
!
        if (impr .ge. 2) then
            ifm = iunifi('MESSAGE')
            call jeimpo(ifm, '&&OP0151.DOMMAGE', 'DOMMAGE')
        end if
!
! ---------------------------------------------------------------------
! ---- FATIGUE VIBRATOIRE
! ---------------------------------------------------------------------
!
    else if (typcal(1:13) .eq. 'FATIGUE_VIBR') then
! ---   STRUCTURE RESULTAT CONTENANT LES CHAM_ELEMS DE SIGMA EQUIVALENT
        call getvid('HISTOIRE', 'RESULTAT', iocc=1, scal=nomres, nbret=nval)
        nomsd = nomres
        call getvid('HISTOIRE', 'MODE_MECA', iocc=1, scal=nommod, nbret=nval)
        nomsd2 = nommod
!
        if (mdomag .eq. 'WOHLER') then
            if (nomopt(11:14) .ne. 'SIGM') then
                call utmess('F', 'FATIGUE1_29', sk=nomopt)
            end if
            pheno = 'FATIGUE'
            call rccome(nommat, pheno, icodre)
            if (icodre .eq. 1) then
                call utmess('F', 'FATIGUE1_88', sk='WOHLER')
            end if
            call rcpare(nommat, pheno, 'WOHLER', icodre)
            if (icodre .eq. 1) then
                call utmess('F', 'FATIGUE1_88', sk='WOHLER')
            end if
        end if
!
! --- CONTRAINTE STATIQUE
!
        call jelira(nomsd//'.ORDR', 'LONUTI', nbord)
        if (nbord .gt. 1) then
            call utmess('F', 'FATIGUE1_84', si=nbord)
        end if
!
        call jenonu(jexnom(nomsd//'.DESC', nomsym), numsym)
        if (numsym .eq. 0) then
            valk(1) = nomsym
            valk(2) = nomsd
            call utmess('F', 'PREPOST4_5', nk=2, valk=valk)
        end if
        call jeveuo(jexnum(nomsd//'.TACH', numsym), 'L', ivch)
        chelrs = zk24(ivch) (1:19)
        if (chelrs .eq. ' ') then
            valk(1) = chelrs
            valk(2) = nomsym
            valk(3) = nomsd
            call utmess('F', 'PREPOST4_6', nk=3, valk=valk)
        end if
        call jeveuo(chelrs//'.CELK', 'L', jcelk)
        ligrel = zk24(jcelk-1+1) (1:19)
        call jelira(chelrs//'.CELV', 'LONMAX', nval)
!
!  -      IL Y A NTCMP COMPOSANTES DANS L OPTION XXXX_EQUI_YYYY
        nbpt = nval/ntcmp
!
! --- CONTRAINTE MODALE
!
        call jenonu(jexnom(nomsd2//'.DESC', nomsym), numsym)
        if (numsym .eq. 0) then
            valk(1) = nomsym
            valk(2) = nomsd2
            call utmess('F', 'PREPOST4_5', nk=2, valk=valk)
        end if
        call jeveuo(jexnum(nomsd2//'.TACH', numsym), 'L', ivch)
        chelrs = zk24(ivch) (1:19)
        if (chelrs .eq. ' ') then
            valk(1) = chelrs
            valk(2) = nomsym
            valk(3) = nomsd2
            call utmess('F', 'PREPOST4_6', nk=3, valk=valk)
        end if
        call jeveuo(chelrs//'.CELK', 'L', jcelk)
        ligrel = zk24(jcelk-1+1) (1:19)
        call jelira(chelrs//'.CELV', 'LONMAX', nval)
!
!  -      IL Y A NTCMP COMPOSANTES DANS L OPTION XXXX_EQUI_YYYY
        nbpt2 = nval/ntcmp
!
        if (nbpt .ne. nbpt2) then
            vali(1) = nbpt
            vali(2) = nbpt2
            call utmess('F', 'FATIGUE1_85', ni=2, vali=vali)
        end if
!
!-- NOMBRE ET NUMERO D ORDRE
        call getvis('HISTOIRE', 'NUME_MODE', iocc=1, nbval=0, nbret=nbord)
        call getvr8('HISTOIRE', 'FACT_PARTICI', iocc=1, nbval=0, nbret=nbord2)
!
        if (nbord .ne. nbord2) then
            call utmess('F', 'FATIGUE1_86')
        end if
!
        nbord = -nbord
        call wkvect('&&OP0151.LMODE', 'V V I', nbord, jordr)
        call getvis('HISTOIRE', 'NUME_MODE', iocc=1, nbval=nbord, vect=zi(jordr), &
                    nbret=ibid)
        call wkvect('&&OP0151.CMODE', 'V V R', nbord, jcoef)
        call getvr8('HISTOIRE', 'FACT_PARTICI', iocc=1, nbval=nbord, vect=zr(jcoef), &
                    nbret=ibid)
!
        if (impr .ge. 1) then
            typeq = nomgde
            typdg = 'CONTRAINTE'
            if (nomopt(6:9) .eq. 'ELNO') typoi = 'NOEUDS'
            if (nomopt(6:9) .eq. 'ELGA') typoi = 'POINTS DE GAUSS'
            valk(1) = typdg
            valk(2) = typoi
            valk(3) = typeq
            vali(1) = nbpt
            vali(2) = nbord
            call utmess('I', 'FATIGUE1_81', nk=3, valk=valk, ni=2, &
                        vali=vali)
        end if
!
        call wkvect('&&OP0151.DOMMAGE', 'V V R', nbpt, ivdmg)
        call dmgmod(nomsym, nomsd, nomsd2, nommat, nbord, &
                    jordr, jcoef, nbpt, ntcmp, numcmp, &
                    impr, zr(ivdmg))
!
        call jedetr('&&OP0151.LMODE')
        call jedetr('&&OP0151.CMODE')
!
    end if
!
! ----- TRANSFORMATION DU VECTEUR DOMMAGE EN UN VRAI CHAM_ELEM
!       ----------------------------------------------------------------
!       ON ALLOUE LE CHAM_ELEM AVEC LA ROUTINE ALCHML
!       PUIS ON RECOPIE DANS .CELV LES VALEURS CALCULEES
!
    chelem = nomu
!
    option = 'TOU_INI_'//nomopt(6:9)
    call alchml(ligrel, option, 'PDOMMAG', 'G', chelem, &
                iret, ' ')
    ASSERT(iret .eq. 0)
!
!
    call jeveuo(chelem//'.CELV', 'E', vr=celv)
    call jelira(chelem//'.CELV', 'LONMAX', ibid)
    ASSERT(ibid .eq. nbpt)
    do i = 1, nbpt
        celv(i) = zr(ivdmg+i-1)
    end do
    call jedetr('&&OP0151.DOMMAGE')
!
7777 continue
!
    call jedema()
end subroutine
