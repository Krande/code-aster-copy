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

subroutine paqmai(nomsd, nomu, nommai, nommet, nomcri, &
                  nomfor, grdvie, forvie, forcri, fordef, &
                  typcha, proaxe, instic, inscri, prec)
! person_in_charge: van-xuan.tran at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/loisem.h"
#include "asterc/lor8em.h"
#include "asterc/r8prem.h"
#include "asterfort/alchml.h"
#include "asterfort/anacri.h"
#include "asterfort/avgrma.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cesexi.h"
#include "asterfort/cesred.h"
#include "asterfort/deltau.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedisp.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jerazo.h"
#include "asterfort/jeveuo.h"
#include "asterfort/reliem.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/utmess.h"
#include "asterfort/vampli.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
    character(len=8) :: nomsd, nomu, nommai, inscri
    character(len=16) :: nommet, nomcri, typcha, proaxe, nomfor, forvie
    character(len=16) :: forcri, grdvie
    aster_logical :: fordef
    real(kind=8) :: instic, prec
! ---------------------------------------------------------------------
! BUT: DETERMINER LE PLUS PETIT CERCLE CIRCONSCRIT AUX POINTS
!      REPRESANTANT LE VECTEUR DE CISAILLEMENT TAU DANS LE PLAN u, v.
! ---------------------------------------------------------------------
! ARGUMENTS:
! NOMSD      IN    K8 : NOM DE LA STRUCTURE DE DONNEES RESULTAT.
! NOMU       IN    K8 : NOM UTILISATEUR DU CALCUL EN FATIGUE.
! NOMMAI     IN    K8 : NOM UTILISATEUR DU MAILLAGE.
! NOMMET     IN    K16: NOM DE LA METHODE DE CALCUL DU CERCLE
!                       CIRCONSCRIT.
! NOMCRI     IN    K16: NOM DU CRITERE.
! NOMFOR   IN    K16: LE NOM DE FORMULE DE GRNADUER EQUIVALENTE
! GRDVIE   IN    K16: NOM DE LA COURBE GRANDEUR EQUIVALENT _DUREE VIE
! TYPCHA     IN    K16: TYPE DE CHARGEMENT (PERIODIQUE OU NON).
! PROAXE     IN    K16: TYPE DE PROJECTION (UN OU DEUX AXES).
!-----------------------------------------------------------------------
    integer(kind=8) :: ibid, ierd, lordr, nbordr, ndim, iret
    integer(kind=8) :: nbma, nbpgt, nbpgmx, jnbpg, ima, tdisp(1), jrwork, tpaq
    integer(kind=8) :: nbpaq, numpaq, nmapaq, nbcmp, bormax, nbpmax
    integer(kind=8) :: nmaini, nbmap, tspaq, iordr, jad, tord(1)
    integer(kind=8) :: jsigv, jsigd, jsigl, imap, nbpg, ipg, icmp, iret1
    integer(kind=8) :: jepsv, jepsd, jepsl, paract(35), jepped, jeppel
    integer(kind=8) :: jepspv, jepspd, jepspl, iret2, jeppev, valep
    integer(kind=8) :: i, kwork, sompgw, sompgs, sompgi, jmail, jgrma
    integer(kind=8) :: n, ninit, nbpggm, nbmagm, nmemo, nncp
    integer(kind=8) :: vali(2), decal, ordini, k, jinst, iord
!
    real(kind=8) :: r8b, val1
!
    complex(kind=8) :: c16b
!
    character(len=1) :: kbid
    character(len=4) :: lsig(6), leps(6)
    character(len=8) :: k8b, motcle(1), tymocl(1)
    character(len=16) :: typres, nomopt
    character(len=19) :: lismai
    character(len=19) :: cesr, ligre, celbid, chsig, chsigs, ces1, ces2
    character(len=19) :: cheps, ces3, ces4, cheppe
    character(len=19) :: chepsp, ces5, ces6, ces7, ces8
    aster_logical :: lbid, crsigm, crepst, crepse, crepsp, creppe
    integer(kind=8), pointer :: nume_ordre(:) => null()
    integer(kind=8), pointer :: paqma(:) => null()
    integer(kind=8), pointer :: cesd(:) => null()
!
!-----------------------------------------------------------------------
!234567                                                              012
!-----------------------------------------------------------------------
    data lsig/'SIXX', 'SIYY', 'SIZZ', 'SIXY', 'SIXZ', 'SIYZ'/
!
    data leps/'EPXX', 'EPYY', 'EPZZ', 'EPXY', 'EPXZ', 'EPYZ'/
!-----------------------------------------------------------------------
!
    call jemarq()
!
! RECUPERATION DU TYPE DE CALCUL MECANIQUE EFFECTUE
!
    call dismoi('TYPE_RESU', nomsd, 'RESULTAT', repk=typres)
    if ((typres(1:9) .ne. 'EVOL_ELAS') .and. (typres(1:9) .ne. 'EVOL_NOLI')) then
        call utmess('F', 'PREPOST4_26')
    end if
!
! CONSTRUCTION DU CHAMP SIMPLE DESTINE A RECEVOIR LES RESULTATS :
! DTAUM,....
!
    call rsexch('F', nomsd, 'SIEF_ELGA', 1, chsig, &
                iret)
!
    call dismoi('NOM_LIGREL', chsig, 'CHAM_ELEM', repk=ligre)
    cesr = '&&PAQMAI.FACY'
    celbid = '&&PAQMAI.BID'
    call alchml(ligre, 'TOU_INI_ELGA', 'PFACY_R', 'V', celbid, &
                ierd, ' ')
    if (ierd .ne. 0) then
        call utmess('A', 'FATIGUE1_1')
    end if
    call celces(celbid, 'V', cesr)
!
! RECUPERATION DU NOMBRE DE NUMEROS D'ORDRE ET DE LA LISTE
! DES NUMEROS D'ORDRE
!
    call rsorac(nomsd, 'TOUT_ORDRE', ibid, r8b, k8b, &
                c16b, r8b, k8b, tord, 1, &
                nbordr)
    lordr = tord(1)
!
    if (nbordr .lt. 0) then
        ndim = -nbordr
    else if (nbordr .gt. 0) then
        ndim = nbordr
    end if
!
    AS_ALLOCATE(vi=nume_ordre, size=ndim)
    call rsorac(nomsd, 'TOUT_ORDRE', ibid, r8b, k8b, &
                c16b, r8b, k8b, nume_ordre, ndim, &
                nbordr)
    if (nume_ordre(1) .eq. 0) then
        call utmess('A', 'PREPOST4_27')
        nbordr = nbordr-1
    end if
!
    ordini = 1
    do k = 2, nbordr
        iord = nume_ordre(k)
        call rsadpa(nomsd, 'L', 1, 'INST', iord, &
                    0, sjv=jinst, styp=kbid)
        if (instic .gt. r8prem()) then
            if (inscri .eq. 'ABSOLU') then
                if (abs(zr(jinst)-instic) .lt. prec) then
                    ordini = k
                    goto 410
                end if
            else
                if (inscri .eq. 'RELATIF') then
                    if (abs(zr(jinst)/instic-1.d0) .lt. prec) then
                        ordini = k
                        goto 410
                    end if
                end if
            end if
        end if
    end do
410 continue
!
    if ((ordini .eq. 1) .and. ((inscri .eq. 'ABSOLU') .or. (inscri .eq. 'RELATIF'))) then
        call utmess('A', 'PREPOST4_48')
    end if
!
!  INITIALISER
    crsigm = .false.
    crepst = .false.
    crepse = .false.
    crepsp = .false.
!---    ANALYSER LE CRITERE
    call anacri(nomcri, nomfor, typcha, 'NON', paract, &
                lbid, crsigm, crepst, crepse, crepsp)
!
!
! IF CRITERE CONTIENT DEFORMATION ELASTIQUE
    creppe = .false.
    if (crepse) then
        if (.not. crepst) then
            call utmess('A', 'PREPOST4_45')
            crepst = .true.
        end if
        if ((.not. crepsp)) then
            call utmess('A', 'PREPOST4_46')
            creppe = .true.
        end if
!
    end if
! RECUPERATION DU NOMBRE DE MAILLES ET DU NOMBRE DE POINTS DE GAUSS
! PAR MAILLE
!
! CAS OU L'ON CALCULE LA FATIGUE SUR TOUT LE MAILLAGE
    call rsexch('F', nomsd, 'SIEF_ELGA', 1, chsig, &
                iret)
    chsigs = '&&PAQMAI.SIELGA'
    call celces(chsig, 'V', chsigs)
    call jeveuo(chsigs(1:19)//'.CESD', 'L', vi=cesd)
    nbma = cesd(1)
!
! CAS OU L'ON CALCULE LA FATIGUE SUR UN OU DES GROUPES DE MAILLES
!
    call wkvect('&&PAQMAI.NBMAGR', 'V V I', 1, jgrma)
    call jerazo('&&PAQMAI.NBMAGR', 1, 1)
    nbmagm = 0
!
    if (nommai .ne. '        ') then
        lismai = '&&PAQMAI.L_MAILLES'
        motcle(1) = 'GROUP_MA'
        tymocl(1) = 'GROUP_MA'
        call reliem(' ', nommai, 'NU_MAILLE', ' ', 0, &
                    1, motcle, tymocl, lismai, nbmagm)
        call jeveuo(lismai, 'L', jmail)
!
!        VECTEUR CONTENANT LE NBR. DE PT. DE GAUSS DES MAILLES DU OU
!        DES GROUPE(S) DE MAILLES
        call jedetr('&&PAQMAI.NBMAGR')
        call wkvect('&&PAQMAI.NBMAGR', 'V V I', nbmagm, jgrma)
!
        do i = 1, nbmagm
            zi(jgrma-1+i) = zi(jmail-1+i)
        end do
    end if
!
!     VECTEUR CONTENANT LE NBR. DE PT. DE GAUSS DES MAILLES DU MAILLAGE
    call wkvect('&&PAQMAI.NBPG', 'V V I', nbma, jnbpg)
!
! NBPGMX : NOMBRE DE POINTS DE GAUSS DANS LES ELEMENTS
!          QUI EN ONT LE PLUS, (EX : ELEMENT 3D = 27)
! NBPGT  : NOMBRE TOTAL DE POINTS DE GAUSS DANS LE MAILLAGE
! NBPGGM : NOMBRE DE POINTS DE GAUSS DANS LE GROUPE DE MAILLES
!
    nbpgmx = 0
    nbpgt = 0
    nbpggm = 0
!
    do ima = 1, nbma
        zi(jnbpg-1+ima) = cesd(5+4*(ima-1)+1)
        nbpgt = nbpgt+cesd(5+4*(ima-1)+1)
        if (cesd(5+4*(ima-1)+1) .gt. nbpgmx) then
            nbpgmx = cesd(5+4*(ima-1)+1)
        end if
    end do
    if (nommai .ne. '        ') then
        do ima = 1, nbmagm
            nbpggm = nbpggm+cesd(5+4*(zi(jmail+ima-1)-1)+ &
                                 1)
        end do
        write (6, *) 'NOMBRE DE POINTS DE GAUSS DU GROUPE DE MAILLES ==>',&
     &              nbpggm
        write (6, *) ' '
    end if
    write (6, *) 'NOMBRE TOTAL DE POINTS DE GAUSS A EXPLORER ==>', nbpgt
    write (6, *) ' '
!
    write (6, *) 'NUMERO DU PAQUET DE MAILLES  -  '//&
     &           'NOMBRE DE POINTS DE GAUSS TRAITES'
!
! CONSTRUCTION DES PAQUETS DE MAILLES.
!
! 1/ DIMENSIONNEMENT DU VECTEUR DE TRAVAIL (RWORK) ET DU VECTEUR
!    CONTENANT LES CARACTERISTIQUES DES PAQUETS DE MAILLES (PAQMA).
!    JEDISP REND LA DIMENSION EN ENTIERS, ON LA CONVERTIT A L'AIDE
!    DES FONCTIONS ENVIMA POUR ALLOUER UN TABLEAU DE REELS.
    call jedisp(1, tdisp)
    tdisp(1) = (tdisp(1)/lor8em())*loisem()
    tdisp(1) = int(0.6d0*tdisp(1))
    call wkvect('&&PAQMAI.RWORK', 'V V R', tdisp(1), jrwork)
!
!       IF (( NOMCRI(1:16) .EQ. 'FATESOCI_MODI_AV' ) .OR.
!      &     FORDEF )THEN
!          NBCMP = 12
!       ELSE
!          NBCMP = 6
!       ENDIF
    nbcmp = 18
!     POUR CALCULER LE NB DE PAQUETS MAXI NOUS NE FAISONS PAS LA MEME
!     CHOSE QUE DANS paqnoe.f PARCE QUE BORMAX EST NATURELLEMENT
!     SURDIMENSIONNEE CAR NBMA TIENT COMPTE DES MAILLES NON VOLUMIQUES.
    bormax = nbma*nbpgmx*nbordr*nbcmp
    val1 = dble(tdisp(1))/dble(bormax)
!
    if (val1 .lt. 1.0d0) then
        nbpmax = int(1.0d0/val1)+1
    else
        nbpmax = 2
    end if
    AS_ALLOCATE(vi=paqma, size=nbpmax*4)
!
! TPAQ   = TAILLE DU PAQUET DE MAILLES
! NBPAQ  = NOMBRE DE PAQUET(S) DE MAILLES
! NUMPAQ = NUMERO DU PAQUET DE MAILLES
! NMAPAQ = NOMBRE DE MAILLES DANS LE PAQUET DE MAILLES
! ZI(JNBPAQ + (NUMPAQ-1)*4 + 2) = NUMERO DE LA MAILLE INITIALE
!                                 DANS LE PAQUET
!
    tpaq = 0
    nbpaq = 0
    numpaq = 0
    nmapaq = 0
!
    do ima = 1, nbma
        tpaq = tpaq+zi(jnbpg-1+ima)*nbordr*nbcmp
        nmapaq = nmapaq+1
!
        if (tpaq .lt. tdisp(1)) then
            if (ima .eq. nbma) then
                numpaq = numpaq+1
                paqma(1+(numpaq-1)*4) = numpaq
                paqma(1+(numpaq-1)*4+1) = tpaq
                paqma(1+(numpaq-1)*4+2) = ima-(nmapaq-1)
                paqma(1+(numpaq-1)*4+3) = nmapaq
                nbpaq = numpaq
            end if
!
        else if ((tpaq .ge. tdisp(1)) .and. (ima .lt. 3)) then
            vali(1) = tdisp(1)
            vali(2) = tpaq
            call utmess('F', 'PREPOST5_67', ni=2, vali=vali)
!
! 2/ STOCKAGE DES NUMEROS DES PAQUETS, DE LA TAILLE DES PAQUETS,
!    DU NUMERO DE LA PREMIERE MAILLE DE CHAQUE PAQUET DE MAILLES,
!    DU NOMBRE DE MAILLE DE CHAQUE PAQUET ET DU NOMBRE DE PAQUET.
!
        else if ((tpaq .ge. tdisp(1)) .and. (ima .gt. 2)) then
! ON RECULE DE DEUX MAILLES POUR ETRE SUR DE NE PAS DEBORDER DU VECTEUR
! DE TRAVAIL (JRWORK).
!
            tpaq = tpaq-zi(jnbpg-1+ima)*nbordr*nbcmp
            tpaq = tpaq-zi(jnbpg-1+ima-1)*nbordr*nbcmp
!
            numpaq = numpaq+1
            paqma(1+(numpaq-1)*4) = numpaq
            paqma(1+(numpaq-1)*4+1) = tpaq
            paqma(1+(numpaq-1)*4+2) = (ima-nmapaq+1)
            paqma(1+(numpaq-1)*4+3) = nmapaq-2
            nbpaq = numpaq
!
            tpaq = zi(jnbpg-1+ima-1)*nbordr*nbcmp
            tpaq = tpaq+zi(jnbpg-1+ima)*nbordr*nbcmp
            nmapaq = 2
            if (ima .eq. nbma) then
                numpaq = numpaq+1
                paqma(1+(numpaq-1)*4) = numpaq
                paqma(1+(numpaq-1)*4+1) = tpaq
                paqma(1+(numpaq-1)*4+2) = ima-(nmapaq-1)
                paqma(1+(numpaq-1)*4+3) = nmapaq
                nbpaq = numpaq
            end if
        end if
!
    end do
!
    if (nbpaq .gt. nbpmax) then
        vali(1) = nbpmax
        vali(2) = nbpaq
        call utmess('F', 'PREPOST5_68', ni=2, vali=vali)
    end if
!
! TRAITEMENT DES PAQUETS DE MAILLES.
!
!  <<REMPLISSAGE>> DU VECTEUR DE TRAVAIL
!
    sompgi = 0
    sompgs = 0
    nmemo = 0
!
    do numpaq = 1, nbpaq
        call jerazo('&&PAQMAI.RWORK', tdisp(1), 1)
        tpaq = paqma(1+(numpaq-1)*4+1)
        nmaini = paqma(1+(numpaq-1)*4+2)
        nbmap = paqma(1+(numpaq-1)*4+3)
        tspaq = tpaq/nbordr
!
!        PERMET D'INITIALISER SOMPGS A CHAQUE PAQUET
        if (numpaq .gt. 1) then
            sompgi = sompgs
        end if
!
        do iordr = 1, nbordr
            if ((numpaq .gt. 1) .and. (iordr .eq. 1)) then
                ninit = nmemo
            else if ((numpaq .eq. 1) .and. (iordr .eq. 1)) then
                ninit = nmaini
            end if
            n = ninit
!
! IF CRITERE CONTIENT CONTRAINTE
            if (crsigm) then
!
                call rsexch('F', nomsd, 'SIEF_ELGA', iordr, chsig, &
                            iret)
!
                ces1 = '&&PAQMAI.SIG_S1'
                ces2 = '&&PAQMAI.SIG_ORDO'
                call celces(chsig, 'V', ces1)
                call cesred(ces1, 0, [ibid], 6, lsig, &
                            'V', ces2)
                call jeexin(ces2(1:19)//'.CESV', iret)
                if (iret .eq. 0) then
                    call utmess('F', 'PREPOST4_29')
                end if
                call jeveuo(ces2(1:19)//'.CESD', 'L', jsigd)
                call jeveuo(ces2(1:19)//'.CESL', 'L', jsigl)
                call jeveuo(ces2(1:19)//'.CESV', 'L', jsigv)
            end if
!
! IF CRITERE CONTIENT DEFORMATION TOTALE
            if (crepst) then
!
                call rsexch('F', nomsd, 'EPSI_ELGA', iordr, cheps, &
                            iret1)
!
                ces3 = '&&PAQMAI.EPS_S3'
                ces4 = '&&PAQMAI.EPS_ORDO'
                call celces(cheps, 'V', ces3)
                call cesred(ces3, 0, [ibid], 6, leps, &
                            'V', ces4)
                call jeexin(ces4(1:19)//'.CESV', iret)
                if (iret .eq. 0) then
                    call utmess('F', 'PREPOST4_34')
                end if
                call jeveuo(ces4(1:19)//'.CESD', 'L', jepsd)
                call jeveuo(ces4(1:19)//'.CESL', 'L', jepsl)
                call jeveuo(ces4(1:19)//'.CESV', 'L', jepsv)
            end if
!
! IF CRITERE CONTIENT DEFORMATION PLASTIQUE
            if (crepsp) then
!
                call rsexch('F', nomsd, 'EPSP_ELGA', iordr, chepsp, &
                            iret2)
!
                ces5 = '&&PAQMAI.EPSP_S3'
                ces6 = '&&PAQMAI.EPSP_ORDO'
                call celces(chepsp, 'V', ces5)
                call cesred(ces5, 0, [ibid], 6, leps, &
                            'V', ces6)
                call jeexin(ces5(1:19)//'.CESV', iret)
                if (iret .eq. 0) then
                    call utmess('F', 'PREPOST4_37')
                end if
                call jeveuo(ces6(1:19)//'.CESD', 'L', jepspd)
                call jeveuo(ces6(1:19)//'.CESL', 'L', jepspl)
                call jeveuo(ces6(1:19)//'.CESV', 'L', jepspv)
            end if
!
! IF CRITERE CONTIENT DEFORMATION ELASTIQUE
            if (creppe) then
!
                call rsexch(' ', nomsd, 'EPSP_ELGA', iordr, cheppe, &
                            valep)
                if (valep .ne. 0) then
                    call utmess('A', 'PREPOST4_46')
                end if
                if (valep .eq. 0) then
                    ces7 = '&&PAQMAI.EPSPE_S3'
                    ces8 = '&&PAQMAI.EPSPE_ORDO'
                    call celces(cheppe, 'V', ces7)
                    call cesred(ces7, 0, [ibid], 6, leps, &
                                'V', ces8)
                    call jeexin(ces7(1:19)//'.CESV', iret)
                    if (iret .eq. 0) then
                        call utmess('F', 'PREPOST4_37')
                    end if
                    call jeveuo(ces8(1:19)//'.CESD', 'L', jepped)
                    call jeveuo(ces8(1:19)//'.CESL', 'L', jeppel)
                    call jeveuo(ces8(1:19)//'.CESV', 'L', jeppev)
                end if
            end if
!
!
            if (numpaq .eq. 1) then
                sompgs = 0
            else if (numpaq .gt. 1) then
                sompgs = sompgi
            end if
            sompgw = 0
            kwork = 0
            decal = 18
!
            do imap = nmaini, nmaini+(nbmap-1)
                if ((imap .gt. nmaini) .and. (numpaq .eq. 1)) then
                    sompgs = sompgs+zi(jnbpg+imap-2)
                    kwork = 1
                    sompgw = sompgw+zi(jnbpg+imap-2)
                end if
!
                if ((imap .gt. nmaini) .and. (numpaq .gt. 1)) then
                    kwork = 1
                    sompgw = sompgw+zi(jnbpg+imap-2)
                end if
!
                if (numpaq .gt. 1) then
                    sompgs = sompgs+zi(jnbpg+imap-2)
                end if
                nbpg = zi(jnbpg+imap-1)
!
                if ((nommai .ne. '        ') .and. (imap .ne. zi(jgrma+n-1))) then
                    n = n-1
                else
                    do ipg = 1, nbpg
!
! BOUCLE SUR LES CONTRAINTES (6 COMPOSANTES)
                        if (crsigm) then
!
                            do icmp = 1, 6
                                call cesexi('C', jsigd, jsigl, imap, ipg, &
                                            1, icmp, jad)
                                if (jad .le. 0) then
                                    if (icmp .eq. 5) then
                                        call utmess('F', 'FATIGUE1_2', si=icmp)
                                    else
                                        call utmess('F', 'PREPOST4_30')
                                    end if
                                else
                                    zr(jrwork+(icmp-1)+(ipg-1)* &
                                       decal+kwork*sompgw*decal+( &
                                       iordr-1)*tspaq) = zr(jsigv-1+ &
                                                            jad)
                                end if
                            end do
                        end if
!
! BOUCLE SUR LES DEFORMATIONS TOTALES (6 COMPOSANTES)
                        if (crepst) then
                            do icmp = 1, 6
                                call cesexi('C', jepsd, jepsl, imap, ipg, &
                                            1, icmp, jad)
                                if (jad .le. 0) then
                                    if (icmp .eq. 5) then
                                        call utmess('F', 'FATIGUE1_3', si=icmp)
                                    else
                                        call utmess('F', 'PREPOST4_35')
                                    end if
                                else
                                    zr(jrwork+(icmp+6-1)+(ipg-1)* &
                                       decal+kwork*sompgw*decal+( &
                                       iordr-1)*tspaq) = zr(jepsv-1+ &
                                                            jad)
                                end if
                            end do
                        end if
!
! BOUCLE SUR LES DEFORMATIONS PLASTIQUES (6 COMPOSANTES)
                        if (crepsp) then
                            do icmp = 1, 6
                                call cesexi('C', jepspd, jepspl, imap, ipg, &
                                            1, icmp, jad)
                                if (jad .le. 0) then
                                    if (icmp .eq. 5) then
                                        call utmess('F', 'FATIGUE1_3', si=icmp)
                                    else
                                        call utmess('F', 'PREPOST4_35')
                                    end if
                                else
                                    zr(jrwork+(icmp+6+6-1)+(ipg- &
                                                            1)*decal+kwork*sompgw*decal+( &
                                       iordr-1)*tspaq) = zr(jepspv-1+ &
                                                            jad)
                                end if
                            end do
                        end if
!
! BOUCLE SUR LES DEFORMATIONS PLASTIQUES (6 COMPOSANTES)
                        if (crepse) then
                            if (valep .eq. 0) then
                                do icmp = 1, 6
                                    call cesexi('C', jepped, jeppel, imap, ipg, &
                                                1, icmp, jad)
                                    if (jad .le. 0) then
                                        if (icmp .eq. 5) then
                                            call utmess('F', 'FATIGUE1_3', si=icmp)
                                        else
                                            call utmess('F', 'PREPOST4_35')
                                        end if
                                    else
                                        zr(jrwork+(icmp+6+6-1)+(ipg- &
                                                                1)*decal+kwork*sompgw*decal+ &
                                           (iordr-1)*tspaq) = zr( &
                                                           jeppev-1+jad)
                                    end if
                                end do
                            else
                                do icmp = 1, 6
                                    zr(jrwork+(icmp+6+6-1)+(ipg-1) &
                                       *decal+kwork*sompgw*decal+( &
                                       iordr-1)*tspaq) = 0.d0
                                end do
                            end if
                        end if
!
                    end do
                end if
                if ((nommai .ne. '        ') .and. (n .lt. nbmagm)) then
                    n = n+1
                end if
!
            end do
            nmemo = n
        end do
!
!         ENDIF
!
        if (nomcri(1:11) .eq. 'VMIS_TRESCA') then
            nomopt = 'DOMA_ELGA'
            call vampli(zr(jrwork), tdisp(1), zi(jnbpg), nbpgt, nbordr, &
                        nmaini, nbmap, tspaq, nomopt, cesr)
            goto 200
        end if
!
        if (typcha .eq. 'PERIODIQUE') then
            call deltau(jrwork, jnbpg, nbpgt, nbordr, ordini, &
                        nmaini, nbmap, numpaq, tspaq, nommet, &
                        nomcri, nomfor, grdvie, forvie, forcri, &
                        cesr)
!
        else if (typcha .eq. 'NON_PERIODIQUE') then
            call avgrma(zr(jrwork), tdisp(1), zi(jnbpg), nbpgt, nbordr, &
                        nmaini, nbmap, numpaq, tspaq, nomcri, &
                        nomfor, grdvie, forvie, fordef, proaxe, &
                        cesr)
        end if
!
200     continue
    end do
!
!
! TRANSFORMATION D'UN CHAM_ELEM SIMPLE EN CHAM_ELEM
!
    call rsexch('F', nomsd, 'SIEF_ELGA', 1, chsig, &
                iret)
    call dismoi('NOM_LIGREL', chsig, 'CHAM_ELEM', repk=ligre)
    call cescel(cesr, ligre, 'TOU_INI_ELGA', ' ', 'NON', &
                nncp, 'G', nomu, 'F', ibid)
!
! MENAGE
!
    call detrsd('CHAM_ELEM', celbid)
    call detrsd('CHAM_ELEM_S', cesr)
    call detrsd('CHAM_ELEM_S', chsigs)
!
    if (crsigm) then
        call detrsd('CHAM_ELEM_S', ces1)
        call detrsd('CHAM_ELEM_S', ces2)
    end if
!
    if (crepst) then
        call detrsd('CHAM_ELEM_S', ces3)
        call detrsd('CHAM_ELEM_S', ces4)
    end if
!
    if (crepsp) then
        call detrsd('CHAM_ELEM_S', ces5)
        call detrsd('CHAM_ELEM_S', ces6)
    end if
!
    if (((creppe)) .and. (valep .eq. 0)) then
        call detrsd('CHAM_ELEM_S', ces7)
        call detrsd('CHAM_ELEM_S', ces8)
    end if
!
    AS_DEALLOCATE(vi=nume_ordre)
    call jedetr('&&PAQMAI.NBMAGR')
    call jedetr('&&PAQMAI.NBPG')
    call jedetr('&&PAQMAI.RWORK')
    AS_DEALLOCATE(vi=paqma)
    if (nommai .ne. '        ') then
        call jedetr('&&PAQMAI.L_MAILLES')
    end if
!
    call jedema()
end subroutine
