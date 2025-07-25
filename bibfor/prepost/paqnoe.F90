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

subroutine paqnoe(nomsd, nomu, nommai, nommet, nomcri, &
                  nomfor, grdvie, forvie, forcri, fordef, &
                  typcha, proaxe, instic, inscri, prec)
! person_in_charge: van-xuan.tran at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/loisem.h"
#include "asterc/lor8em.h"
#include "asterc/r8prem.h"
#include "asterfort/anacri.h"
#include "asterfort/avgrno.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnscno.h"
#include "asterfort/cnscre.h"
#include "asterfort/cnsred.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/dtauno.h"
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
! BUT: CONSTRUIRE LES PAQUETS DE NOEUDS AFIN DE CALCULER LE VECTEUR
!      CISAILLEMENT TAU DANS LE PLAN u, v.
! ---------------------------------------------------------------------
! ARGUMENTS:
! NOMSD    IN    K8 : NOM DE LA STRUCTURE DE DONNEES RESULTAT.
! NOMU     IN    K8 : NOM UTILISATEUR DU CALCUL EN FATIGUE.
! NOMMAI   IN    K8 : NOM UTILISATEUR DU MAILLAGE.
! NOMMET   IN    K16: NOM DE LA METHODE DE CALCUL DU CERCLE
!                     CIRCONSCRIT.
! NOMCRI   IN    K16: NOM DU CRITERE D'ENDOMMAGEMENT PAR FATIGUE.
! NOMFOR   IN    K16: LE NOM DE FORMULE DE GRNADUER EQUIVALENTE
! GRDVIE   IN    K16: NOM DE LA COURBE GRANDEUR EQUIVALENT _DUREE VIE
! TYPCHA   IN    K16: TYPE DE CHARGEMENT (PERIODIQUE OU NON).
! NOMCRI   IN    K16: NOM DU CRITERE AVEC PLANS CRITIQUES.
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ibid, lordr, nbordr, ndim, iret, iret1
    integer(kind=8) :: nbno, ino, tdisp(1), jrwork, tpaq, iret2, iret3
    integer(kind=8) :: nbpaq, numpaq, nnopaq, bormax, nbpmax, nbp0, bor0
    integer(kind=8) :: nbcmp, k, ordini, iord, jinst
    integer(kind=8) :: nnoini, nbnop, tspaq, iordr
    integer(kind=8) :: jsigv, jsigd, jsigl, jepsd, jepsl, jepsv, kwork
    integer(kind=8) :: jepspd, jepspl, jepspv, iret4, jepped, jeppel
    integer(kind=8) :: inop, somnow, icmp, jnoeu, nunoe, jeppev
    integer(kind=8) :: vali(2), paract(35), valep, ic, tord(1)
    aster_logical :: post, crsigm, crepst, crepse, crepsp, lbid, creppe
!
    real(kind=8) :: r8b, val1, val2, resu(7), vresu(24), valpar(35)
    complex(kind=8) :: c16b
!
    character(len=1) :: kbid
    character(len=4) :: lsig(6), leps(6)
    character(len=8) :: k8b, lresu(24), motcle(2), tymocl(2)
    character(len=16) :: typres, nomopt
    character(len=19) :: cnsr, lisnoe, cheppe, cns7, cns8
    character(len=19) :: chsig, cheps, cns1, cns2, cns3, cns4, chsig1, chsig2
    character(len=19) :: chepsp, cns5, cns6
    integer(kind=8), pointer :: nume_ordre(:) => null()
    integer(kind=8), pointer :: paqno(:) => null()
!
!
!-----------------------------------------------------------------------
!234567                                                              012
!-----------------------------------------------------------------------
    data lsig/'SIXX', 'SIYY', 'SIZZ', 'SIXY', 'SIXZ', 'SIYZ'/
!
    data leps/'EPXX', 'EPYY', 'EPZZ', 'EPXY', 'EPXZ', 'EPYZ'/
!
    data lresu/'DTAUM1', 'VNM1X', 'VNM1Y', 'VNM1Z', 'SINMAX1',&
     &             'SINMOY1', 'EPNMAX1', 'EPNMOY1', 'SIGEQ1', 'NBRUP1',&
     &             'ENDO1', 'DTAUM2', 'VNM2X', 'VNM2Y', 'VNM2Z',&
     &             'SINMAX2', 'SINMOY2', 'EPNMAX2', 'EPNMOY2', 'SIGEQ2',&
     &             'NBRUP2', 'ENDO2', 'VMIS', 'TRESCA'/
!-----------------------------------------------------------------------
!
    call jemarq()
!
!               1234567890123456789
    chsig = '                   '
    chsig1 = '                   '
    chsig2 = '                   '
    cheps = '                   '
    chepsp = '                   '
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
    cnsr = '&&PAQNOE.FACY'
    call cnscre(nommai, 'FACY_R', 24, lresu, 'V', &
                cnsr)
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
    if (nume_ordre(1) .eq. 0) then
        call utmess('A', 'PREPOST4_27')
        nbordr = nbordr-1
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
! IF CRITERE CONTIENT DEFORMATION ELASTIQUE
    creppe = .false.
    if (crepse) then
        if (.not. crepst) then
            call utmess('A', 'PREPOST4_45')
            crepst = .true.
        end if
!
        if ((.not. crepsp)) then
            call utmess('A', 'PREPOST4_46')
            creppe = .true.
        end if
!
    end if
!
! CREATION D'UN OBJET JEVEUX CONTENANT LA LISTE DES NUMEROS
! DE NOEUDS AINSI QUE LE NOMBRE DE NOEUDS
!
    lisnoe = '&&PAQNOE.L_NOEUDS'
    motcle(1) = 'GROUP_MA'
    tymocl(1) = 'GROUP_MA'
    motcle(2) = 'GROUP_NO'
    tymocl(2) = 'GROUP_NO'
    call reliem(' ', nommai, 'NU_NOEUD', ' ', 0, &
                2, motcle, tymocl, lisnoe, nbno)
    call jeveuo(lisnoe, 'L', jnoeu)
!
    write (6, *) 'NOMBRE TOTAL DE NOEUDS A TRAITER ==>', nbno
    write (6, *) ' '
    write (6, *) 'NUMERO DU PAQUET DE NOEUDS   -   '//&
     &           'NOMBRE DE NOEUDS TRAITES'
!
! CONSTRUCTION DES PAQUETS DE NOEUDS.
!
! 1/ DIMENSIONNEMENT DU VECTEUR DE TRAVAIL (RWORK) ET DU VECTEUR
!    CONTENANT LES CARACTERISTIQUES DES PAQUETS DE NOEUDS (PAQNO).
!    JEDISP REND LA DIMENSION EN ENTIERS, ON LA CONVERTIT A L'AIDE
!    DES FONCTIONS ENVIMA POUR ALLOUER UN TABLEAU DE REELS
    call jedisp(1, tdisp)
    tdisp(1) = (tdisp(1)*loisem())/lor8em()
    tdisp(1) = int(0.6d0*tdisp(1))
    call wkvect('&&PAQNOE.RWORK', 'V V R', tdisp(1), jrwork)
!
    nbcmp = 18
!
    bormax = nbno*nbordr*nbcmp
    val1 = dble(tdisp(1))/dble(bormax)
!
!     ON TIENT COMPTE DU RECUL DE DEUX NOEUDS POUR CALCULER NBPMAX
    if (val1 .lt. 1.0d0) then
        nbp0 = int(1.0d0/val1)+1
        bor0 = nbp0*2*nbordr*nbcmp
        val2 = dble(tdisp(1))/dble(bormax+bor0)
        nbpmax = int(1.0d0/val2)+1
    else
        nbpmax = 2
    end if
    AS_ALLOCATE(vi=paqno, size=nbpmax*4)
!
! TPAQ   = TAILLE DU PAQUET DE NOEUDS
! NBPAQ  = NOMBRE DE PAQUET(S) DE NOEUDS
! NUMPAQ = NUMERO DU PAQUET DE NOEUDS
! NNOPAQ = NOMBRE DE NOEUDS DANS LE PAQUET DE NOEUDS
! ZI(JPAQNO + (NUMPAQ-1)*4 + 2) = NUMERO DU NOEUD INITIAL DANS LE PAQUET
!
    tpaq = 0
    nbpaq = 0
    numpaq = 0
    nnopaq = 0
!
    do ino = 1, nbno
        tpaq = tpaq+nbordr*nbcmp
        nnopaq = nnopaq+1
!
        if (tpaq .lt. tdisp(1)) then
            if (ino .eq. nbno) then
                numpaq = numpaq+1
                paqno(1+(numpaq-1)*4) = numpaq
                paqno(1+(numpaq-1)*4+1) = tpaq
                paqno(1+(numpaq-1)*4+2) = ino-(nnopaq-1)
                paqno(1+(numpaq-1)*4+3) = nnopaq
                nbpaq = numpaq
            end if
!
        else if ((tpaq .ge. tdisp(1)) .and. (ino .lt. 3)) then
            vali(1) = tdisp(1)
            vali(2) = tpaq
            call utmess('F', 'PREPOST5_67', ni=2, vali=vali)
!
! 2/ STOCKAGE DES NUMEROS DES PAQUETS, DE LA TAILLE DES PAQUETS,
!    DU NUMERO DE LA PREMIERE MAILLE DE CHAQUE PAQUET DE MAILLES,
!    DU NOMBRE DE MAILLE DE CHAQUE PAQUET ET DU NOMBRE DE PAQUET.
!
        else if ((tpaq .ge. tdisp(1)) .and. (ino .gt. 2)) then
! ON RECULE DE DEUX NOEUDS POUR ETRE SUR DE NE PAS DEBORDER DU VECTEUR
! DE TRAVAIL (JRWORK).
!
            tpaq = tpaq-(2*nbordr*nbcmp)
            numpaq = numpaq+1
            paqno(1+(numpaq-1)*4) = numpaq
            paqno(1+(numpaq-1)*4+1) = tpaq
            paqno(1+(numpaq-1)*4+2) = ino-(nnopaq-1)
            paqno(1+(numpaq-1)*4+3) = nnopaq-2
            nbpaq = numpaq
!
            tpaq = 2*nbordr*nbcmp
            nnopaq = 2
            if (ino .eq. nbno) then
                numpaq = numpaq+1
                paqno(1+(numpaq-1)*4) = numpaq
                paqno(1+(numpaq-1)*4+1) = tpaq
                paqno(1+(numpaq-1)*4+2) = ino-(nnopaq-1)
                paqno(1+(numpaq-1)*4+3) = nnopaq
                nbpaq = numpaq
            end if
        end if
!
    end do
!
    if (nbpaq .gt. nbpmax) then
        vali(1) = nbpmax
        vali(2) = nbpaq
        call utmess('F', 'PREPOST5_70', ni=2, vali=vali)
    end if
!
! TRAITEMENT DES PAQUETS DE NOEUDS.
!
!  <<REMPLISSAGE>> DU VECTEUR DE TRAVAIL
!
    do numpaq = 1, nbpaq
        call jerazo('&&PAQNOE.RWORK', tdisp(1), 1)
        tpaq = paqno(1+(numpaq-1)*4+1)
        nnoini = paqno(1+(numpaq-1)*4+2)
        nbnop = paqno(1+(numpaq-1)*4+3)
        tspaq = tpaq/nbordr
!
        do iordr = 1, nbordr
!
! IF CRITERE CONTIENT CONTRAINTE
            if (crsigm) then
!
                call rsexch(' ', nomsd, 'SIGM_NOEU', iordr, chsig1, &
                            iret1)
                call rsexch(' ', nomsd, 'SIEF_NOEU', iordr, chsig2, &
                            iret2)
!
                if ((iret1 .ne. 0) .and. (iret2 .ne. 0)) then
                    call utmess('A', 'PREPOST4_38')
                end if
!
                if (iret1 .eq. 0) then
                    chsig = chsig1
                else if (iret2 .eq. 0) then
                    chsig = chsig2
                end if
!
                cns1 = '&&PAQNOE.SIG_S1'
                cns2 = '&&PAQNOE.SIG_ORDO'
                call cnocns(chsig, 'V', cns1)
                call cnsred(cns1, 0, [ibid], 6, lsig, &
                            'V', cns2)
                call jeexin(cns2(1:19)//'.CNSV', iret)
                if (iret .eq. 0) then
                    call utmess('F', 'PREPOST4_40')
                end if
                call jeveuo(cns2(1:19)//'.CNSD', 'L', jsigd)
                call jeveuo(cns2(1:19)//'.CNSL', 'L', jsigl)
                call jeveuo(cns2(1:19)//'.CNSV', 'L', jsigv)
            end if
!
! IF CRITERE CONTIENT DEFORMATION TOTALE
            if (crepst) then
!
                call rsexch('F', nomsd, 'EPSI_NOEU', iordr, cheps, &
                            iret3)
!
                cns3 = '&&PAQNOE.EPS_S3'
                cns4 = '&&PAQNOE.EPS_ORDO'
                call cnocns(cheps, 'V', cns3)
                call cnsred(cns3, 0, [ibid], 6, leps, &
                            'V', cns4)
                call jeexin(cns4(1:19)//'.CNSV', iret)
                if (iret .eq. 0) then
                    call utmess('F', 'PREPOST4_41')
                end if
                call jeveuo(cns4(1:19)//'.CNSD', 'L', jepsd)
                call jeveuo(cns4(1:19)//'.CNSL', 'L', jepsl)
                call jeveuo(cns4(1:19)//'.CNSV', 'L', jepsv)
!
            end if
!
! IF CRITERE CONTIENT DEFORMATION PLASTQIUE
            if (crepsp) then
!
                call rsexch('F', nomsd, 'EPSP_NOEU', iordr, chepsp, &
                            iret4)
!
                cns5 = '&&PAQNOE.EPSP_S5'
                cns6 = '&&PAQNOE.EPSP_ORDO'
                call cnocns(chepsp, 'V', cns5)
                call cnsred(cns5, 0, [ibid], 6, leps, &
                            'V', cns6)
                call jeexin(cns6(1:19)//'.CNSV', iret)
                if (iret .eq. 0) then
                    call utmess('F', 'PREPOST4_43')
                end if
                call jeveuo(cns6(1:19)//'.CNSD', 'L', jepspd)
                call jeveuo(cns6(1:19)//'.CNSL', 'L', jepspl)
                call jeveuo(cns6(1:19)//'.CNSV', 'L', jepspv)
!
            end if
!
! IF CRITERE CONTIENT DEFORMATION ELASTIQUE
            if (creppe) then
!
                call rsexch(' ', nomsd, 'EPSP_NOEU', iordr, cheppe, &
                            valep)
                if (valep .ne. 0) then
                    call utmess('A', 'PREPOST4_46')
                end if
                if (valep .eq. 0) then
                    cns7 = '&&PAQNOE.EPSPE_S7'
                    cns8 = '&&PAQNOE.EPSPE_ORDO'
                    call cnocns(cheppe, 'V', cns7)
                    call cnsred(cns7, 0, [ibid], 6, leps, &
                                'V', cns8)
                    call jeexin(cns8(1:19)//'.CNSV', iret)
                    if (iret .eq. 0) then
                        call utmess('F', 'PREPOST4_43')
                    end if
                    call jeveuo(cns8(1:19)//'.CNSD', 'L', jepped)
                    call jeveuo(cns8(1:19)//'.CNSL', 'L', jeppel)
                    call jeveuo(cns8(1:19)//'.CNSV', 'L', jeppev)
                end if
            end if
!
            kwork = 0
            somnow = 0
!
            do inop = nnoini, nnoini+(nbnop-1)
                if (inop .gt. nnoini) then
                    kwork = 1
                    somnow = somnow+1
                end if
!
                nunoe = zi(jnoeu+inop-1)
!
! BOUCLE SUR LES CONTRAINTES (6 COMPOSANTES)
                if (crsigm) then
                    do icmp = 1, 6
                        if (zl(jsigl+(icmp-1)+(nunoe-1)*6)) then
                            zr(jrwork+(icmp-1)+kwork*somnow*18+ &
                               (iordr-1)*tspaq) = zr(jsigv+(icmp-1)+ &
                                                     (nunoe-1)*6)
                        else if (icmp .eq. 5) then
                            call utmess('F', 'FATIGUE1_2', si=icmp)
                        else
                            call utmess('F', 'PREPOST4_30')
                        end if
                    end do
                end if
!
! BOUCLE SUR LES DEFORMATIONS TOTALES (6 COMPOSANTES)
                if (crepst) then
                    do icmp = 1, 6
                        if (zl(jepsl+(icmp-1)+(nunoe-1)*6)) then
                            zr(jrwork+(icmp+6-1)+kwork*somnow*18 &
                               +(iordr-1)*tspaq) = zr(jepsv+(icmp-1) &
                                                      +(nunoe-1)*6)
                        else if (icmp .eq. 5) then
                            call utmess('F', 'FATIGUE1_3', si=icmp)
                        else
                            call utmess('F', 'PREPOST4_35')
                        end if
                    end do
                end if
!
! BOUCLE SUR LES DEFORMATIONS TOTALES (6 COMPOSANTES)
                if (crepsp) then
                    do icmp = 1, 6
                        if (zl(jepspl+(icmp-1)+(nunoe-1)*6)) then
                            zr(jrwork+(icmp+6+6-1)+kwork*somnow* &
                               18+(iordr-1)*tspaq) = zr(jepspv+( &
                                                        icmp-1)+(nunoe-1)*6)
                        else if (icmp .eq. 5) then
                            call utmess('F', 'FATIGUE1_3', si=icmp)
                        else
                            call utmess('F', 'PREPOST4_35')
                        end if
                    end do
                end if
!
                if (creppe) then
                    if (valep .eq. 0) then
!
                        do icmp = 1, 6
                            if (zl(jeppel+(icmp-1)+(nunoe-1)*6)) then
                                zr(jrwork+(icmp+6+6-1)+kwork* &
                                   somnow*18+(iordr-1)*tspaq) = &
                                    zr(jeppev+(icmp-1)+(nunoe-1)*6)
                            else if (icmp .eq. 5) then
                                call utmess('F', 'FATIGUE1_3', si=icmp)
                            else
                                call utmess('F', 'PREPOST4_35')
                            end if
                        end do
!
                    else
                        do ic = 1, 6
                            zr(jrwork+(ic+6+6-1)+kwork*somnow*18 &
                               +(iordr-1)*tspaq) = 0.d0
                        end do
                    end if
                end if
!
            end do
        end do
!
!
!
!
        if (nomcri(1:11) .eq. 'VMIS_TRESCA') then
            nomopt = 'DOMA_NOEUD'
            call vampli(zr(jrwork), tdisp(1), zi(jnoeu), nbno, nbordr, &
                        nnoini, nbnop, tspaq, nomopt, cnsr)
            goto 200
        end if
!
        if (typcha .eq. 'PERIODIQUE') then
            post = .false.
!
            call dtauno(jrwork, zi(jnoeu), nbno, nbordr, ordini, &
                        nnoini, nbnop, tspaq, nommet, nomcri, &
                        nomfor, grdvie, forvie, forcri, nommai, &
                        cnsr, k8b, post, valpar, vresu)
!
        else if (typcha .eq. 'NON_PERIODIQUE') then
!
!   POUR POST_FATIGUE
            post = .false.
!
            call avgrno(zr(jrwork), tdisp(1), zi(jnoeu), nbno, nbordr, &
                        nnoini, nbnop, tspaq, nomcri, nomfor, &
                        grdvie, forvie, fordef, nommai, proaxe, &
                        k8b, cnsr, post, resu)
        end if
!
200     continue
    end do
!
! TRANSFORMATION D'UN CHAM_NO SIMPLE EN CHAM_NO
!
    call cnscno(cnsr, ' ', 'NON', 'G', nomu, &
                'F', ibid)
!
! MENAGE
    call detrsd('CHAM_NO_S', cnsr)
!
    if (crsigm) then
        call detrsd('CHAM_NO_S', cns1)
        call detrsd('CHAM_NO_S', cns2)
    end if
!
    if (crepst) then
        call detrsd('CHAM_NO_S', cns3)
        call detrsd('CHAM_NO_S', cns4)
    end if
!
    if (crepsp) then
        call detrsd('CHAM_NO_S', cns5)
        call detrsd('CHAM_NO_S', cns6)
    end if
!
    if (creppe .and. (valep .eq. 0)) then
        call detrsd('CHAM_NO_S', cns7)
        call detrsd('CHAM_NO_S', cns8)
    end if
!
    AS_DEALLOCATE(vi=nume_ordre)
    call jedetr('&&PAQNOE.RWORK')
    AS_DEALLOCATE(vi=paqno)
!
    call jedema()
end subroutine
