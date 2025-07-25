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
subroutine meca01(optio0, nbordr, jordr, nchar, jcha, &
                  kcha, tbgrca, resuco, resuc1, leres1, &
                  noma, modele, ligrmo, mate, cara, &
                  chvarc, codret)
! person_in_charge: josselin.delmas at edf.fr
! ----------------------------------------------------------------------
! COMMANDE DE CALC_ERREUR SPECIFIQUE AUX INDICATEURS D'ERREUR
! ----------------------------------------------------------------------
! IN  OPTIO0 : OPTION A TRAITER
! IN  NBORDR : NOMBRE DE NUMEROS D'ORDRE
! IN  JORDR  : ADRESSES DES NUMEROS D'ORDRES
! IN  NCHAR  : NOMBRE DE CHARGES
! IN  JCHA   : ADRESSES DES CHARGES
! IN  KCHA   : NOM JEVEUX OU SONT STOCKEES LES CHARGES
! IN  TBGRCA : TABLEAU DES GRANDEURS CARACTERISTIQUES (HM)
! IN  RESUCO : NOM DE CONCEPT RESULTAT
! IN  RESUC1 : NOM DE CONCEPT DE LA COMMANDE CALC_ERREUR
! IN  LERES1 : NOM DE CONCEPT RESULTAT A ENRICHIR
! IN  NOMA   : NOM DU MAILLAGE
! IN  MODELE : NOM DU MODELE
! IN  LIGRMO : LISTE DES GROUPES DU MODELE
! IN  MATE   : NOM DU CHAMP MATERIAU
! IN  CARA   : NOM DU CHAMP DES CARACTERISTIQUES ELEMENTAIRES
! OUT CODRET : CODE DE RETOUR AVEC 0 SI TOUT VA BIEN
!              1 : PROBLEMES DE DONNEES
!              2 : PROBLEMES DE RESULTATS
! ----------------------------------------------------------------------
!
    implicit none
!
!     --- ARGUMENTS ---
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/alchml.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/erglob.h"
#include "asterfort/ernozz.h"
#include "asterfort/exisd.h"
#include "asterfort/exithm.h"
#include "asterfort/exixfe.h"
#include "asterfort/gettco.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jerecu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/ltnotb.h"
#include "asterfort/lxlgut.h"
#include "asterfort/mecact.h"
#include "asterfort/mecara.h"
#include "asterfort/mechti.h"
#include "asterfort/medom1.h"
#include "asterfort/qintzz.h"
#include "asterfort/qires1.h"
#include "asterfort/reslgn.h"
#include "asterfort/reslo2.h"
#include "asterfort/resloc.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexc1.h"
#include "asterfort/rsexc2.h"
#include "asterfort/rsnoch.h"
#include "asterfort/tbliva.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
!
    integer(kind=8) :: nbordr, jordr, nchar, jcha
    integer(kind=8) :: codret
    real(kind=8) :: tbgrca(3)
    character(len=8) :: noma, resuco, resuc1, modele, cara
    character(len=19) :: kcha, chvarc
    character(len=19) :: leres1
    character(len=24) :: ligrmo
    character(len=24) :: mate
    character(len=*) :: optio0
!
!
!     --- VARIABLES LOCALES ---
!
    character(len=6) :: nompro
    parameter(nompro='MECA01')
!
    integer(kind=8) :: npacri
    parameter(npacri=2)
!
    integer(kind=8) :: iordr, jfin, jaux, tabido(5)
    integer(kind=8) :: np, nd, ncharp, nchard, jchap, jchad
    integer(kind=8) :: iret, iret1, iad
    integer(kind=8) :: iainst, ii
    integer(kind=8) :: iaux, ibid
    integer(kind=8) :: vali
    integer(kind=8) :: irxfem
!
    real(kind=8) :: rbid, rundf, theta, deltat
    real(kind=8) :: time, erp, erd, s, longc, presc
!
    character(len=2) :: cret
    character(len=8) :: k8b
    character(len=8) :: ctype
    character(len=8) :: resup, resud
    character(len=16) :: option, optiop, optiod, lipacr(npacri), tysd
    character(len=19) :: kchap, kchad, tabp, tabd
    character(len=24) :: blan24, k24b
    character(len=24) :: chs
    character(len=24) :: chdepp, chsgpn, chsgdn, mateco
    character(len=24) :: chsig, chsigp, chsigd, chsign
    character(len=24) :: chsigm, chdepm, cherrm, chsigx
    character(len=24) :: chcara(18), chelem, chtime
    character(len=24) :: cherre, cherrn
    character(len=24) :: ligrch, ligrcp, ligrcd
    character(len=24) :: chvois, cvoisx
    character(len=24) :: valk(2)
!
    complex(kind=8) :: cbid, valc
    aster_logical :: yaxfem, yathm, perman
!
!=======================================================================
! 1. PREALABLES
!=======================================================================
!
    rbid = 0.d0
    cbid = (0.d0, 0.d0)
    codret = 0
    rundf = r8vide()
!
! 1.1. ==> VERIFICATIONS
!
! 1.1.1. ==> L'OPTION
!
    iaux = lxlgut(optio0)
    if (iaux .gt. 16) then
        call utmess('F', 'INDICATEUR_98', sk=optio0(1:iaux))
    else
        option = '                '
        option(1:iaux) = optio0(1:iaux)
    end if
!
    do ii = 1, nchar
        ligrch = zk8(jcha+ii-1)//'.CHME.LIGRE'
    end do
!
! 1.1.2. ==> LE TYPE DE SD
!
    call gettco(resuco, tysd)
!
! 1.2. ==> INITIALISATIONS
!
    blan24 = ' '
!               12   345678   9012345678901234
    kchap = '&&'//nompro//'.CHARGESP  '
    kchad = '&&'//nompro//'.CHARGESD  '
    chsigp = blan24
    chsigd = blan24
    chsigx = blan24
    ncharp = 0
    nchard = 0
    yaxfem = .false.
!
!=======================================================================
! 2. OPTION "ERME_ELEM"
!=======================================================================
!
    if (option .eq. 'ERME_ELEM') then
!
! 2.1. ==> PREALABLES
!
!--- RECHERCHE DES VOISINS
!--- (CHGEOM RECHERCHE A PARTIR DU MODELE ET PAS DES CHARGES)
        call reslo2(modele, ligrmo, chvois, cvoisx, tabido)
! --- EST-CE DE LA THM ?
        call exithm(modele, yathm, perman)
! --- EST-CE DU XFEM ?
        call exixfe(modele, irxfem)
        if (irxfem .ne. 0) yaxfem = .true.
!
! --- POUR DE LA THM EN TRANSITOIRE, ON DEVRA RECUPERER LES INFORMATIONS
!     DU PAS DE TEMPS PRECEDENT
        jfin = 1
        if (yathm) then
            if (.not. perman) then
                jfin = 2
            end if
        end if
!
!--- INITIALISATIONS  : DEPLACEMENTS ET CONTRAINTES
!
        chdepm = ' '
        chdepp = ' '
        chsigm = ' '
        chsigp = ' '
!
! 2.2. ==> BOUCLE SUR LES NUMEROS D'ORDRE
!
        do iaux = 1, nbordr
!
            call jemarq()
            call jerecu('V')
            iordr = zi(jordr+iaux-1)
!
! 2.2.1 ==> SAISIT ET VERIFIE LA COHERENCE DES DONNEES MECANIQUES
!           RECUPERE LES CHARGES POUR LE NUMERO D'ORDRE IORDR
            call medom1(modele, mate, mateco, cara, kcha, &
                        nchar, resuco, iordr)
            call jeveuo(kcha//'.LCHA', 'L', jcha)
            call mecara(cara, chcara)
!
            if ((tysd .eq. 'EVOL_ELAS') .or. (tysd .eq. 'EVOL_NOLI')) then
!--- RECUPERATION DES INSTANTS CORRESPONDANT A IORDR ET IORDR-1
                do jaux = 1, jfin
!
                    ibid = iordr+1-jaux
!
                    if (ibid .lt. 0) then
                        call utmess('I', 'INDICATEUR_3', sk=nompro)
                        goto 299
                    end if
!
                    call rsadpa(resuco, 'L', 1, 'INST', ibid, &
                                0, sjv=iainst, styp=k8b)
!
                    if (jaux .eq. 1) then
                        time = zr(iainst)
                        deltat = rundf
                        theta = rundf
                    else
                        deltat = time-zr(iainst)
!
! - --RECUPERATION DU PARM_THETA CORRESPONDANT A IORDR
                        call jenonu(jexnom(resuco//'           .NOVA', 'PARM_THETA'), iad)
                        if (iad .eq. 0) then
                            theta = 0.57d0
                            call utmess('A', 'INDICATEUR_4', sk=resuco)
                        else
                            call rsadpa(resuco, 'L', 1, 'PARM_THETA', iordr, &
                                        0, sjv=iad, styp=k8b)
                            theta = zr(iad)
                            if ((theta .gt. 1.d0) .or. (theta .lt. 0.d0)) then
                                call utmess('F', 'INDICATEUR_5', sk=resuco)
                            end if
                        end if
! - --
                    end if
                end do
            else
                time = 0.d0
                deltat = rundf
                theta = rundf
            end if
!--- CREATION DE LA CARTE DES INSTANTS
            call mechti(noma, time, deltat, theta, chtime)
!
! 2.2.2 ==> RECUPERATION DES CHAMPS DE CONTRAINTES AUX NOEUDS
            do jaux = 1, jfin
! VERIFIE L'EXISTENCE DU CHAMP
! S'IL EXISTE ON RECUPERE SON NOM SYMBOLIQUE
                ibid = iordr+1-jaux
                call rsexc2(1, 2, resuco, 'SIGM_ELNO', ibid, &
                            k24b, option, iret)
                call rsexc2(2, 2, resuco, 'SIEF_ELNO', ibid, &
                            k24b, option, iret)
!
!--- SI AUCUN CHAMP N'EXISTE, ON SORT
                if (iret .gt. 0) goto 299
!
! 2.2.3. ==> VERIFIE SI LE CHAMP EST CALCULE SUR TOUT LE MODELE
                call dismoi('NOM_LIGREL', k24b, 'CHAM_ELEM', repk=ligrch)
                if (ligrch .ne. ligrmo) then
                    call codent(ibid, 'G', k8b)
                    valk(1) = option
                    valk(2) = k8b
                    call utmess('A', 'INDICATEUR_2', nk=2, valk=valk)
                    goto 299
                end if
!--- ARCHIVAGE DU NOM DU CHAMP DE CONTRAINTES
                if (jaux .eq. 1) then
                    chsigp = k24b
                else
                    chsigm = k24b
                end if
            end do
!
!--- RECUPERATION DU CHAMP DE CONTRAINTES AUX NOEUDS PAR SOUS
!    ELEMENT (OPTION 'SIEF_SENO_SEGA'), CAS X-FEM UNIQUEMENT
            if (yaxfem) then
                call rsexc2(1, 1, resuco, 'SISE_ELNO', iordr, &
                            chsigx, option, iret)
                if (iret .gt. 0) goto 299
            end if
! ---------------------------------------------------------------------
! 2.2.4 ==> POUR DE LA THM : ON RECUPERE ...
! ---------------------------------------------------------------------
            if (yathm) then
! -----------------------------
! 2.2.4.1. LES CHAMPS DE DEPLACEMENTS
! -----------------------------
                do jaux = 1, jfin
!
                    ibid = iordr+1-jaux
                    call rsexc2(1, 1, resuco, 'DEPL', ibid, &
                                k24b, option, iret1)
                    if (iret1 .gt. 0) then
                        call codent(ibid, 'G', k8b)
                        valk(1) = resuco
                        valk(2) = k8b
                        call utmess('A', 'CALCULEL3_11', nk=2, valk=valk)
                        goto 299
                    end if
                    if (jaux .eq. 1) then
                        chdepp = k24b
                    else
                        chdepm = k24b
                    end if
                end do
! ---------------------------------------
! 2.2.4.2. LES GRANDEURS CARACTERISTIQUES
! ---------------------------------------
!
                longc = tbgrca(1)
                presc = tbgrca(2)
                if (longc .le. 0.d0 .or. presc .le. 0.d0) then
                    call utmess('F', 'INDICATEUR_28')
                end if
! -----------------------------
! 2.2.4.1. LE CHAMP D'ESTIMATEURS A L'INSTANT PRECEDENT
! -----------------------------
!
                if (.not. perman) then
!
                    if (iordr .eq. 1) then
!
! INSTANT INITIAL : CREATION D'UN CHAM_ELEM NUL
!                         12   345678   9012345678901234
                        cherrm = '&&'//nompro//'_ERREUR_M       '
                        call alchml(ligrmo, option, 'PERREM', 'V', cherrm, &
                                    iret, ' ')
                        if (iret .ne. 0) then
                            call utmess('A', 'CALCULEL5_4')
                            goto 299
                        end if
!
                    else
! SINON, ON RECUPERE LE CHAMP DE L'INSTANT PRECEDENT
                        ibid = iordr-1
                        call rsexc2(1, 1, resuco, 'ERME_ELEM', ibid, &
                                    k24b, option, iret1)
                        if (iret1 .gt. 0) then
                            call codent(ibid, 'G', k8b)
                            valk(1) = resuco
                            valk(2) = k8b
                            call utmess('F', 'INDICATEUR_24', nk=2, valk=valk)
                            goto 299
                        end if
                        cherrm = k24b
                    end if
!
                end if
!
            end if
!
! 2.2.7. ==> RECUPERE LE NOM SYMBOLIQUE DU CHAMP DE L'OPTION CALCULEE
!            POUR LE NUMERO D'ORDRE IORDR
            call rsexc1(leres1, option, iordr, chelem)
!
! 2.2.8. ==> CALCULE L'ESTIMATEUR D'ERREUR EN RESIDU LOCAL
!
            call resloc(modele, ligrmo, yaxfem, yathm, tbgrca, &
                        perman, chtime, mateco, chsigm, chsigp, &
                        chsigx, chdepm, chdepp, cherrm, zk8(jcha), &
                        nchar, tabido, chvois, cvoisx, chelem)
!
! 2.2.9. ==> VERIFIE L'EXISTENCE DU CHAMP CHELEM
            call exisd('CHAMP_GD', chelem, iret)
!
!--- SI LE CHAMP N'EXISTE PAS, ON SORT
            if (iret .eq. 0) then
                codret = 1
                call jedema()
                goto 999
            end if
!
! 2.2.10. ==> CALCUL DE L'ESTIMATEUR GLOBAL A PARTIR DES ESTIMATEURS
!             LOCAUX
            call erglob(chelem, yathm, perman, option, iordr, &
                        resuco, leres1)
!
! 2.2.11. ==> NOTE LE NOM D'UN CHAMP19 DANS UNE SD_RESULTAT
            call rsnoch(leres1, option, iordr)
!
299         continue
!
            call jedema()
!
        end do
!
!=======================================================================
! 3. OPTION "ERME_ELNO"
!=======================================================================
!
    else if (option .eq. 'ERME_ELNO') then
!
        do 11, iaux = 1, nbordr
!
            call jemarq()
!
            call jerecu('V')
            iordr = zi(jordr+iaux-1)
            if (iordr .eq. 0) then
                call jedema()
                goto 11
            end if
!
            call rsexc2(1, 1, resuco, 'ERME_ELEM', iordr, &
                        cherre, option, iret1)
!
            if (iret1 .eq. 0) then
                call rsexc1(leres1, option, iordr, cherrn)
                call reslgn(ligrmo, option, cherre, cherrn)
                call rsnoch(leres1, option, iordr)
            end if
!
            call jedema()
!
11          continue
!
!=======================================================================
! 4. OPTION "QIRE_ELEM"
!=======================================================================
!
            else if (option .eq. 'QIRE_ELEM') then
!
! 4.1. ==> PREALABLES
! 4.1.1. ==> RECUPERE LES NOMS DES SD RESULTAT
            call getvid(' ', 'RESULTAT', scal=resup, nbret=np)
            call getvid(' ', 'RESU_DUAL', scal=resud, nbret=nd)
!
! 4.1.2. ==> RECHERCHE DES VOISINS
            call reslo2(modele, ligrmo, chvois, cvoisx, tabido)
! 4.1.3. ==>  RECUPERE LES NOMS SYMBOLIQUES DES TABLES
            tabp = ' '
            tabd = ' '
            call ltnotb(resup, 'ESTI_GLOB', tabp)
            call ltnotb(resud, 'ESTI_GLOB', tabd)
!
! 4.2. ==> BOUCLE SUR LES NUMEROS D'ORDRE
!
            do 12, iaux = 1, nbordr
!
                call jemarq()
                iordr = zi(jordr+iaux-1)
!
! 4.2.1. ==> CALCULE LE COEFFICIENT S
!----- RECUPERE ERRE_ABSO DANS LA TABLE A PARTIR DU NUMERO D'ORDRE
                lipacr(1) = 'NUME_ORDR'
                lipacr(2) = 'OPTION'
!
                call tbliva(tabp, npacri, lipacr, [iordr], [rbid], &
                            [cbid], 'ERME_ELEM', 'EGAL', [0.d0], 'ERRE_ABSO', &
                            ctype, vali, erp, valc, valk(1), &
                            iret)
                call tbliva(tabd, npacri, lipacr, [iordr], [rbid], &
                            [cbid], 'ERME_ELEM', 'EGAL', [0.d0], 'ERRE_ABSO', &
                            ctype, vali, erd, valc, valk(1), &
                            iret)
                s = sqrt(erd/erp)
!----- CREE UNE CARTE CONSTANTE
                chs = '&&OP0069.CH_NEUT_R'
                call mecact('V', chs, 'MODELE', ligrmo, 'NEUT_R', &
                            ncmp=1, nomcmp='X1', sr=s)
!
! 4.2.2. ==> SAISIE ET VERIFIE LA COHERENCE DES DONNEES MECANIQUES
                call medom1(modele, mate, mateco, cara, kchap, &
                            ncharp, resup, iordr)
                call medom1(modele, mate, mateco, cara, kchad, &
                            nchard, resud, iordr)
                call jeveuo(kchap//'.LCHA', 'L', jchap)
                call jeveuo(kchad//'.LCHA', 'L', jchad)
!
! 4.2.3. ==> VERIFIE L'EXISTENCE DU CHAMP DANS LE RESUPRIM
!          S'IL EXISTE RECUPERE SON NOM SYMBOLIQUE
                call rsexc2(1, 2, resup, 'SIGM_ELNO', iordr, &
                            chsigp, option, iret)
                call rsexc2(2, 2, resup, 'SIEF_ELNO', iordr, &
                            chsigp, option, iret)
!
!         SI AUCUN CHAMP N'EXISTE, ON SORT
                if (iret .gt. 0) goto 499
!
! 4.2.4. ==> VERIFIE L'EXISTENCE DU CHAMP DANS LE RESUDUAL
!         S'IL EXISTE RECUPERE SON NOM SYMBOLIQUE
                call rsexc2(1, 2, resud, 'SIGM_ELNO', iordr, &
                            chsigd, option, iret)
                call rsexc2(2, 2, resud, 'SIEF_ELNO', iordr, &
                            chsigd, option, iret)
!
!         SI AUCUN CHAMP N'EXISTE, ON SORT
                if (iret .gt. 0) goto 499
!
! 4.2.5. ==> RECUPERE LE NOM DE L'OPTION CALCULEE POUR CHACUN DES CHAMPS
                call dismoi('NOM_OPTION', chsigp, 'CHAM_ELEM', repk=optiop)
                call dismoi('NOM_OPTION', chsigd, 'CHAM_ELEM', repk=optiod)
!
! 4.2.6. ==> VERIFIE SI LE CHAMP EST CALCULE SUR TOUT LE MODELE
                call dismoi('NOM_LIGREL', chsigp, 'CHAM_ELEM', repk=ligrcp)
                call dismoi('NOM_LIGREL', chsigd, 'CHAM_ELEM', repk=ligrcd)
                if (ligrcp .ne. ligrmo .or. ligrcd .ne. ligrmo) then
                    call codent(iordr, 'G', k8b)
                    valk(1) = option
                    valk(2) = k8b
                    call utmess('A', 'INDICATEUR_2', nk=2, valk=valk)
                    goto 499
                end if
!
! 4.2.7. ==> RECUPERE L'ADRESSE JEVEUX DE L'INSTANT DE CALCUL
!          POUR LE NUMERO D'ORDRE IORDR
                if ((tysd .eq. 'EVOL_ELAS') .or. (tysd .eq. 'EVOL_NOLI')) then
                    call rsadpa(resuco, 'L', 1, 'INST', iordr, &
                                0, sjv=iainst, styp=k8b)
                    time = zr(iainst)
                else
                    time = 0.d0
                end if
!
! 4.2.8. ==> CREE UNE CARTE D'INSTANTS
                call mechti(noma, time, rundf, rundf, chtime)
!
! 4.2.9. ==> RECUPERE LE NOM SYMBOLIQUE DU CHAMP DE L'OPTION CALCULEE
!           POUR LE NUMERO D'ORDRE IORDR
                call rsexc1(leres1, option, iordr, chelem)
!
! 4.2.10. ==> CALCULE L'ESTIMATEUR D'ERREUR EN RESIDU LOCAL
!
                call qires1(modele, ligrmo, chtime, chsigp, chsigd, &
                            zk8(jchap), zk8(jchad), ncharp, nchard, chs, &
                            mateco, chvois, tabido, chelem)
!
! 4.2.11. ==> VERIFIE L'EXISTENCE DU CHAMP CHELEM
                call exisd('CHAMP_GD', chelem, iret)
!
!--- SI LE CHAMP N'EXISTE PAS, ON SORT
                if (iret .eq. 0) then
                    codret = 1
                    call jedema()
                    goto 999
                end if
!
! 4.2.12. ==> CALCUL DE L'ESTIMATEUR GLOBAL A PARTIR DES ESTIMATEURS
!             LOCAUX
                call erglob(chelem, .false._1, .false._1, option, iordr, &
                            resuco, leres1)
!
! 4.2.13. ==> NOTE LE NOM D'UN CHAMP19 DANS UNE SD_RESULTAT
                call rsnoch(leres1, option, iordr)
!
499             continue
!
                call jedema()
!
12              continue
!
!=======================================================================
! 5. OPTION "QIRE_ELNO"
!=======================================================================
!
                else if (option .eq. 'QIRE_ELNO') then
!
                do 13, iaux = 1, nbordr
!
                    call jemarq()
!
                    iordr = zi(jordr+iaux-1)
                    call rsexc2(1, 1, resuco, 'QIRE_ELEM', iordr, &
                                cherre, option, iret1)
!
                    if (iret1 .eq. 0) then
                        call rsexc1(leres1, option, iordr, cherrn)
                        call reslgn(ligrmo, option, cherre, cherrn)
                        call rsnoch(leres1, option, iordr)
                    end if
!
                    call jedema()
!
13                  continue
!
!=======================================================================
! 6. OPTIONS "QIZ1_ELEM" ET "QIZ2_ELEM"
!=======================================================================
!
                    else if (option .eq. 'QIZ1_ELEM' .or. option .eq. 'QIZ2_ELEM') then
!
! 6.1. ==> RECUPERE LES NOMS DES SD RESULTAT
                    call getvid(' ', 'RESULTAT', scal=resup, nbret=np)
                    call getvid(' ', 'RESU_DUAL', scal=resud, nbret=nd)
!
! 6.2. ==> BOUCLE SUR LES NUMEROS D'ORDRE
!
                    do 14, iaux = 1, nbordr
!
                        call jemarq()
!
                        iordr = zi(jordr+iaux-1)
!
! 6.2.1. ==> SAISIT ET VERIFIE LA COHERENCE DES DONNEES MECANIQUES
                        call medom1(modele, mate, mateco, cara, kchap, &
                                    ncharp, resup, iordr)
                        call medom1(modele, mate, mateco, cara, kchad, &
                                    nchard, resud, iordr)
                        call jeveuo(kchap//'.LCHA', 'L', jchap)
                        call jeveuo(kchad//'.LCHA', 'L', jchad)
!
! 6.2.2. ==> RECUPERE SON NOM SYMBOLIQUE DU CHAMP DE CONTRAINTES LISSE
!            DANS LE RESUPRIM
                        call rsexc2(1, 1, resup, 'SI'//option(3:4)//'_NOEU', iordr, &
                                    chsgpn, option, iret)
!
! 6.2.3. ==>  RECUPERE SON NOM SYMBOLIQUE DU CHAMP DE CONTRAINTES LISSE
!             DANS LE RESUDUAL
                        call rsexc2(1, 1, resud, 'SI'//option(3:4)//'_NOEU', iordr, &
                                    chsgdn, option, iret)
!
! 6.2.4. ==> RECUPERE SON NOM SYMBOLIQUE DU CHAMP DE CONTRAINTES CALCULE
!            DANS LE RESUPRIM
                        call rsexc2(1, 1, resup, 'SIEF_ELGA', iordr, &
                                    chsigp, option, iret)
                        if (iret .gt. 0) goto 699
!
! 6.2.5 ==> RECUPERE SON NOM SYMBOLIQUE DU CHAMP DE CONTRAINTES CALCULE
!           DANS LE RESUDUAL
                        call rsexc2(1, 1, resud, 'SIEF_ELGA', iordr, &
                                    chsigd, option, iret)
                        if (iret .gt. 0) goto 699
!
! 6.2.6 ==> CALCUL
                        call rsexc1(resuc1, option, iordr, chelem)
!
                        call qintzz(modele, ligrmo, mateco, chsigp, chsigd, &
                                    chsgpn, chsgdn, chelem)
!
! 6.2.7.0. ==> RECUPERE L'ADRESSE JEVEUX DE L'INSTANT DE CALCUL
!          POUR LE NUMERO D'ORDRE IORDR
                        if ((tysd .eq. 'EVOL_ELAS') .or. (tysd .eq. 'EVOL_NOLI')) then
                            call rsadpa(resuco, 'L', 1, 'INST', iordr, &
                                        0, sjv=iainst, styp=k8b)
                            time = zr(iainst)
                        else
                            time = 0.d0
                        end if
!
! 6.2.7.1. ==> RECUPERE LE CHAMP DE VARIABLE DE COMMANDE
                        call vrcins(modele, mate, cara, time, chvarc, &
                                    cret)
!
! 6.2.7.2. ==> CALCUL DE L'ESTIMATEUR GLOBAL A PARTIR DES ESTIMATEURS
!             LOCAUX
                        call ernozz(modele, chsigp, mateco, chsgpn, chvarc, &
                                    option, ligrmo, iordr, resuco, leres1, &
                                    chelem)
!
! 6.2.8. ==> NOTE LE NOM D'UN CHAMP19 DANS UNE SD_RESULTAT
                        call rsnoch(leres1, option, iordr)
!
699                     continue
!
                        call jedema()
!
14                      continue
!
!=======================================================================
! 7. OPTIONS "ERZ1_ELEM" ET "ERZ2_ELEM"
!=======================================================================
!
                        else if (option .eq. 'ERZ1_ELEM' .or. option .eq. 'ERZ2_ELEM') then
!
                        do 15, iaux = 1, nbordr
!
                            call jemarq()
                            call jerecu('V')
!
                            iordr = zi(jordr+iaux-1)
                            call medom1(modele, mate, mateco, cara, kcha, &
                                        nchar, resuco, iordr)
                            call jeveuo(kcha//'.LCHA', 'L', jcha)
                            call mecara(cara, chcara)
                            call rsexc2(1, 1, resuco, 'SIEF_ELGA', iordr, &
                                        chsig, option, iret)
                            if (iret .gt. 0) then
                                call utmess('A', 'CALCULEL3_7', sk=option)
                                codret = 2
                                goto 999
                            end if
!
!
                            call rsexc2(1, 1, resuco, 'SI'//option(3:4)//'_NOEU', iordr, &
                                        chsign, option, iret)
!
                            if (iret .eq. 0) then
                                if ((tysd .eq. 'EVOL_ELAS') .or. (tysd .eq. 'EVOL_NOLI')) then
                                    call rsadpa(resuco, 'L', 1, 'INST', iordr, &
                                                0, sjv=iainst, styp=k8b)
                                    time = zr(iainst)
                                else
                                    time = 0.d0
                                end if
                                call vrcins(modele, mate, cara, time, chvarc, &
                                            cret)
                                call rsexc1(leres1, option, iordr, chelem)
                                call ernozz(modele, chsig, mateco, chsign, chvarc, &
                                            option, ligrmo, iordr, resuco, leres1, &
                                            chelem)
                                call rsnoch(leres1, option, iordr)
                            end if
!
                            call jedema()
!
15                          continue
!
!=======================================================================
! N. OPTION NE CORRESPONDANT PAS AUX INDICATEURS D'ERREUR
!=======================================================================
!
                            else
!
!                 123456   890123456789
                            valk(1) = nompro//'            '
                            valk(2) = option//'   '
                            call utmess('F', 'INDICATEUR_99', nk=2, valk=valk)
!
                            end if
!
999                         continue
!
                            end subroutine
