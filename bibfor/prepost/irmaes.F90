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
!
subroutine irmaes(idfimd, nomaas, nomamd, nbimpr, caimpi, &
                  modnum, nuanom, nomtyp, nnotyp, sdcarm, &
                  carael)
!
    implicit none
!
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/as_mmhcyw.h"
#include "asterfort/as_mmhraw.h"
#include "asterfort/cesexi.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8)  :: nomaas, nomtyp(*), sdcarm, carael
    character(len=64) :: nomamd
    integer(kind=8) :: nbimpr, caimpi(10, nbimpr), modnum(MT_NTYMAX)
    integer(kind=8) :: nnotyp(*), nuanom(MT_NTYMAX, *), nvtyge
    med_idt :: idfimd
!
! --------------------------------------------------------------------------------------------------
!
!  IMPR_RESU - IMPRESSION DANS LE MAILLAGE DES ELEMENTS DE STRUCTURE
!
! --------------------------------------------------------------------------------------------------
!
! IN  :
!   IDFIMD  K*   ENTIER LIE AU FICHIER MED OUVERT
!   NOMAAS  K8   NOM DU MAILLAGE ASTER
!   NOMAMD  K*   NOM DU MAILLAGE MED
!   NBIMPR  K*   NOMBRE D'IMPRESSIONS
!   CAIMPI  K*   ENTIERS POUR CHAQUE IMPRESSION
!   MODNUM  I(*) INDICATEUR SI LA SPECIFICATION DE NUMEROTATION DES
!                NOEUDS DES MAILLES EST DIFFERENTES ENTRE ASTER ET MED
!   NUANOM  I*   TABLEAU DE CORRESPONDANCE DES NOEUDS MED/ASTER.
!                NUANOM(ITYP,J): NUMERO DANS ASTER DU J IEME NOEUD DE LA
!                MAILLE DE TYPE ITYP DANS MED.
!   NOMTYP  K8(*)NOM DES TYPES POUR CHAQUE MAILLE
!   NNOTYP  I(*) NOMBRE DE NOEUDS PAR TYPES DE MAILLES
!   SDCARM  K*   SD_CARA_ELEM EN CHAM_ELEM_S
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: codret, ipoin, ityp, letype, ino, iret, nbcmp, iad, iadzr, iadep
    integer(kind=8) :: ima, nbsect, nbcouc, nbsp, nummai, nbmail, jpoin, ibid, iadex, iadr1
    integer(kind=8) :: jcnxma(MT_NTYMAX), jepama(MT_NTYMAX), joripmf(MT_NTYMAX)
    integer(kind=8) :: jr1eptu(MT_NTYMAX), jorituy(MT_NTYMAX), nmatyp(MT_NTYMAX)
    integer(kind=8) :: j_pmf_cesl, j_pmf_cesd, j_cg_cesl, j_cg_cesd
    integer(kind=8) :: j_tu_cesl, j_tu_cesd, j_tu_or_cesl, j_tu_or_cesd
    integer(kind=8) :: icmp_cq_ep, icmp_cq_ex, icmp_gr_se, icmp_gr_ex, icmp_pf_g1
    integer(kind=8) :: icmp_tu_r1, icmp_tu_ep, icmp_tu_g1, icmp_tu_g2, icmp_tu_g3, icmp_tu_g4
    integer(kind=8) :: nbgamm
    integer(kind=8), parameter :: edfuin = 0, edelst = 5, ednoda = 0
    character(len=3) :: saux03
    character(len=8) :: saux08, nomo
    character(len=64), parameter :: atepai = 'SCALE', atangv = 'ANGLE'
    aster_logical :: exicoq, exituy, exipmf, exigri, okgr, okcq, oktu, okpf
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
    real(kind=8), pointer :: cg_cesv(:) => null()
    real(kind=8), pointer :: tu_cesv(:) => null()
    real(kind=8), pointer :: tu_or_cesv(:) => null()
    real(kind=8), pointer :: pmf_cesv(:) => null()
    character(len=8), pointer :: tmp_cesc(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jeveuo(nomaas//'.TYPMAIL', 'L', vi=typmail)
    call jelira(nomaas//'.TYPMAIL', 'LONMAX', nbmail)
    call jeveuo(jexatr(nomaas//'.CONNEX', 'LONCUM'), 'L', jpoin)
    call jeveuo(nomaas//'.CONNEX', 'L', vi=connex)
!
    exicoq = ASTER_FALSE
    exigri = ASTER_FALSE
    call exisd('CHAM_ELEM_S', sdcarm//'.CARCOQUE', iret)
    if (iret .ne. 0) then
        call jeveuo(sdcarm//'.CARCOQUE  .CESV', 'L', vr=cg_cesv)
        call jeveuo(sdcarm//'.CARCOQUE  .CESL', 'L', j_cg_cesl)
        call jeveuo(sdcarm//'.CARCOQUE  .CESD', 'L', j_cg_cesd)
        call jeveuo(sdcarm//'.CARCOQUE  .CESC', 'L', vk8=tmp_cesc)
!       COQUE  : Si EP     et EXCENT
!       GRILLE : Si SECT_L et DIST_N
        nbcmp = zi(j_cg_cesd+1)
        icmp_cq_ep = indik8(tmp_cesc, 'EP', 1, nbcmp)
        icmp_cq_ex = indik8(tmp_cesc, 'EXCENT', 1, nbcmp)
        if ((icmp_cq_ep .ne. 0) .and. (icmp_cq_ex .ne. 0)) then
            exicoq = ASTER_TRUE
        end if
        icmp_gr_se = indik8(tmp_cesc, 'SECT_L', 1, nbcmp)
        icmp_gr_ex = indik8(tmp_cesc, 'DIST_N', 1, nbcmp)
        if ((icmp_gr_se .ne. 0) .and. (icmp_gr_ex .ne. 0)) then
            exigri = ASTER_TRUE
        end if
    end if
!
    exituy = ASTER_FALSE
    call dismoi('NOM_MODELE', carael, 'CARA_ELEM', ibid, nomo, 'F', iret)
    call dismoi('EXI_TUYAU', nomo, 'MODELE', ibid, saux03, 'F', iret)
    if (saux03 .eq. 'OUI') then
        exituy = ASTER_TRUE
        call jeveuo(sdcarm//'.CARGEOPO  .CESV', 'L', vr=tu_cesv)
        call jeveuo(sdcarm//'.CARGEOPO  .CESL', 'L', j_tu_cesl)
        call jeveuo(sdcarm//'.CARGEOPO  .CESD', 'L', j_tu_cesd)
        call jeveuo(sdcarm//'.CARGEOPO  .CESC', 'L', vk8=tmp_cesc)
        nbcmp = zi(j_tu_cesd+1)
        icmp_tu_r1 = indik8(tmp_cesc, 'R1', 1, nbcmp)
        icmp_tu_ep = indik8(tmp_cesc, 'EP1', 1, nbcmp)
        if ((icmp_tu_r1 .eq. 0) .or. (icmp_tu_ep .eq. 0)) then
            exituy = ASTER_FALSE
        end if
    end if
!
    exipmf = ASTER_FALSE
    call exisd('CHAM_ELEM_S', sdcarm//'.CAFIBR', iret)
    if (iret .ne. 0) then
        exipmf = ASTER_TRUE
    end if
    if (exipmf) then
        call jeveuo(sdcarm//'.CARORIEN  .CESV', 'L', vr=pmf_cesv)
        call jeveuo(sdcarm//'.CARORIEN  .CESL', 'L', j_pmf_cesl)
        call jeveuo(sdcarm//'.CARORIEN  .CESD', 'L', j_pmf_cesd)
        call jeveuo(sdcarm//'.CARORIEN  .CESC', 'L', vk8=tmp_cesc)
        nbcmp = zi(j_pmf_cesd+1)
        icmp_pf_g1 = indik8(tmp_cesc, 'GAMMA', 1, nbcmp)
    end if
!
    nbgamm = 1
    if (exituy) then
        call jeveuo(sdcarm//'.CARORIEN  .CESV', 'L', vr=tu_or_cesv)
        call jeveuo(sdcarm//'.CARORIEN  .CESL', 'L', j_tu_or_cesl)
        call jeveuo(sdcarm//'.CARORIEN  .CESD', 'L', j_tu_or_cesd)
        call jeveuo(sdcarm//'.CARORIEN  .CESC', 'L', vk8=tmp_cesc)
        nbcmp = zi(j_tu_or_cesd+1)
        icmp_tu_g1 = indik8(tmp_cesc, 'GAMMA', 1, nbcmp)
        icmp_tu_g2 = indik8(tmp_cesc, 'GAMMA2', 1, nbcmp)
        icmp_tu_g3 = indik8(tmp_cesc, 'GAMMA3', 1, nbcmp)
        icmp_tu_g4 = indik8(tmp_cesc, 'GAMMA4', 1, nbcmp)
        nbgamm = 4
        if (icmp_tu_g4 .eq. 0) nbgamm = 3
    end if
!
! DECOMPTE DU NOMBRE DE MAILLES PAR TYPE
    nmatyp(:) = 0
    do ima = 1, nbmail
        nmatyp(typmail(ima)) = nmatyp(typmail(ima))+1
    end do
!
!   CREATION D'UN VECTEURS PAR TYPE DE MAILLE PRESENT CONTENANT LA CONNECTIVITE DES MAILLE/TYPE
!   (CONNECTIVITE = NOEUDS + UNE VALEUR BIDON(0) SI BESOIN)
    do ityp = 1, MT_NTYMAX
        if (nmatyp(ityp) .ne. 0) then
            call wkvect('&&IRMAES.CNX.'//nomtyp(ityp), 'V V I', &
                        nnotyp(ityp)*nmatyp(ityp), jcnxma(ityp))
            if (exicoq .or. exigri) then
                ! paramètres : épaisseur, excentrement
                call wkvect('&&IRMAES.EPAI.'//nomtyp(ityp), 'V V R', &
                            2*nmatyp(ityp), jepama(ityp))
            end if
            if (exipmf) then
                ! paramètres : gamma
                call wkvect('&&IRMAES.ORIP.'//nomtyp(ityp), 'V V R', &
                            nmatyp(ityp), joripmf(ityp))
            end if
            if (exituy) then
                ! paramètres : gamma[3 ou 4]
                call wkvect('&&IRMAES.ORIT.'//nomtyp(ityp), 'V V R', &
                            nbgamm*nmatyp(ityp), jorituy(ityp))
                ! paramètres : R1, EP1
                call wkvect('&&IRMAES.RAYO.'//nomtyp(ityp), 'V V R', &
                            2*nmatyp(ityp), jr1eptu(ityp))
            end if
        end if
    end do
!
!   ON PARCOURT TOUTES LES MAILLES. POUR CHACUNE D'ELLES, ON STOCKE SA CONNECTIVITE
!   LA CONNECTIVITE EST FOURNIE EN STOCKANT TOUS LES NOEUDS A LA SUITE POUR UNE MAILLE DONNEE.
!   C'EST CE QU'ON APPELLE LE MODE ENTRELACE DANS MED A LA FIN DE CETTE PHASE, NMATYP CONTIENT
!   LE NOMBRE DE MAILLES POUR CHAQUE TYPE
    nmatyp(:) = 0
    do ima = 1, nbmail
        ityp = typmail(ima)
        ipoin = zi(jpoin-1+ima)
        nmatyp(ityp) = nmatyp(ityp)+1
!
        if (exicoq .or. exigri) then
!           COQUE  : épaisseur,    excentrement
!           GRILLE : épaisseur=0 , excentrement
            iadzr = jepama(ityp)+2*(nmatyp(ityp)-1)
            zr(iadzr) = 0.0d0
            zr(iadzr+1) = 0.0d0
            call cesexi('C', j_cg_cesd, j_cg_cesl, ima, 1, 1, icmp_cq_ep, iadep)
            call cesexi('C', j_cg_cesd, j_cg_cesl, ima, 1, 1, icmp_cq_ex, iadex)
            if ((iadep .gt. 0) .and. (iadex .gt. 0)) then
                zr(iadzr) = cg_cesv(iadep)
                zr(iadzr+1) = cg_cesv(iadex)
            else
                call cesexi('C', j_cg_cesd, j_cg_cesl, ima, 1, 1, icmp_gr_se, iadep)
                call cesexi('C', j_cg_cesd, j_cg_cesl, ima, 1, 1, icmp_gr_ex, iadex)
                if ((iadep .gt. 0) .and. (iadex .gt. 0)) then
                    zr(iadzr) = 0.0d0
                    zr(iadzr+1) = cg_cesv(iadex)
                end if
            end if
        end if
        if (exipmf) then
!           gamma
            iadzr = joripmf(ityp)+nmatyp(ityp)-1
            call cesexi('C', j_pmf_cesd, j_pmf_cesl, ima, 1, 1, icmp_pf_g1, iad)
            if (iad .gt. 0) then
                zr(iadzr) = pmf_cesv(iad)
            else
                zr(iadzr) = 0.0d0
            end if
        end if
        if (exituy) then
            iadzr = jorituy(ityp)+nbgamm*(nmatyp(ityp)-1)
            zr(iadzr:iadzr+2) = 0.0d0
!           gamma1
            call cesexi('C', j_tu_or_cesd, j_tu_or_cesl, ima, 1, 1, icmp_tu_g1, iad)
            if (iad .gt. 0) zr(iadzr) = tu_or_cesv(iad)
!           gamma2
            call cesexi('C', j_tu_or_cesd, j_tu_or_cesl, ima, 1, 1, icmp_tu_g2, iad)
            if (iad .gt. 0) zr(iadzr+1) = tu_or_cesv(iad)
!           gamma3
            call cesexi('C', j_tu_or_cesd, j_tu_or_cesl, ima, 1, 1, icmp_tu_g3, iad)
            if (iad .gt. 0) zr(iadzr+2) = tu_or_cesv(iad)
            if (nbgamm .eq. 4) then
!               gamma4
                call cesexi('C', j_tu_or_cesd, j_tu_or_cesl, ima, 1, 1, icmp_tu_g4, iad)
                if (iad .gt. 0) then
                    zr(iadzr+3) = tu_or_cesv(iad)
                else
                    zr(iadzr+3) = 0.0d0
                end if
            end if
!           R1 (Rextérieur = Rmax), épaisseur
            call cesexi('C', j_tu_cesd, j_tu_cesl, ima, 1, 1, icmp_tu_r1, iadr1)
            call cesexi('C', j_tu_cesd, j_tu_cesl, ima, 1, 1, icmp_tu_ep, iadep)
            iadzr = jr1eptu(ityp)+2*(nmatyp(ityp)-1)
            if ((iadr1 .gt. 0) .and. (iadep .gt. 0)) then
                ! Rmin , Rmax
                zr(iadzr) = tu_cesv(iadr1)-tu_cesv(iadep)
                zr(iadzr+1) = tu_cesv(iadr1)
            else
                zr(iadzr) = 0.d0
                zr(iadzr+1) = 0.d0
            end if
        end if
!       CONNECTIVITE DE LA MAILLE TYPE ITYP DANS VECT CNX
!       I) MAILLES DONT LA NUMEROTATION DES NOEUDS ENTRE ASTER ET MED EST IDENTIQUE.
        if (modnum(ityp) .eq. 0) then
            do ino = 1, nnotyp(ityp)
                ibid = jcnxma(ityp)+(nmatyp(ityp)-1)*nnotyp(ityp)+ino-1
                zi(ibid) = connex(ipoin-1+ino)
            end do
!       II) MAILLES DONT LA NUMEROTATION DES NOEUDS ENTRE ASTER ET MED EST DIFFERENTE (CF LRMTYP).
        else
            do ino = 1, nnotyp(ityp)
                ibid = jcnxma(ityp)+(nmatyp(ityp)-1)*nnotyp(ityp)+ino-1
                zi(ibid) = connex(ipoin-1+nuanom(ityp, ino))
            end do
        end if
!
    end do
!
!   ECRITURE
    cletype: do letype = 1, nbimpr
!       PASSAGE DU NUMERO DE TYPE MED AU NUMERO DE TYPE ASTER
        nbsp = caimpi(3, letype)
        nbcouc = caimpi(4, letype)
        nbsect = caimpi(5, letype)
        nummai = caimpi(6, letype)
        ityp = caimpi(8, letype)
        nvtyge = caimpi(9, letype)
        if (nmatyp(ityp) .eq. 0) cycle cletype
        ! CAS : GRILLE, COQUE, TUYAU, PMF
        okgr = (nummai .eq. 0) .and. (nbcouc .eq. 1) .and. (nbsect .eq. 0) .and. (nbsp .eq. 1)
        okcq = (nummai .eq. 0) .and. (nbcouc .ge. 1) .and. (nbsect .eq. 0) .and. &
               (nbsp .eq. 3*nbcouc)
        oktu = (nummai .eq. 0) .and. (nbcouc .ge. 1) .and. (nbsect .ge. 1)
        okpf = (nummai .ne. 0) .and. (nbcouc .eq. 0) .and. (nbsect .eq. 0)
        if (.not. (okgr .or. okcq .or. oktu .or. okpf)) cycle cletype
!
!       la connectivite est fournie en stockant tous les noeuds a la suite pour une maille.
!       c'est ce que med appelle le mode entrelace
        call as_mmhcyw(idfimd, nomamd, zi(jcnxma(ityp)), nnotyp(ityp)*nmatyp(ityp), edfuin, &
                       nmatyp(ityp), edelst, nvtyge, ednoda, codret)
        if (codret .ne. 0) then
            saux08 = 'mmhcyw'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
!
        if (okcq .or. okgr) then
            ! coque, grille : attribut SCALE ==>  [ epaisseur , excentrement ]
            call as_mmhraw(idfimd, nomamd, nvtyge, atepai, nmatyp(ityp), &
                           zr(jepama(ityp)), codret)
            if (codret .ne. 0) then
                saux08 = 'mmhraw'
                call utmess('F', 'DVP_97', sk=saux08, si=codret)
            end if
        else if (okpf) then
            ! pmf : attribut ANGLE ==> gamma
            call as_mmhraw(idfimd, nomamd, nvtyge, atangv, nmatyp(ityp), &
                           zr(joripmf(ityp)), codret)
            if (codret .ne. 0) then
                saux08 = 'mmhraw'
                call utmess('F', 'DVP_97', sk=saux08, si=codret)
            end if
        else if (oktu) then
            ! tuyau : attribut ANGLE ==> 3 ou 4 gamma
            !       : attribut SCALE ==> [Rmin, Rmax] (calculés à partir de R1 et EP
            call as_mmhraw(idfimd, nomamd, nvtyge, atangv, nmatyp(ityp), &
                           zr(jorituy(ityp)), codret)
            if (codret .ne. 0) then
                saux08 = 'mmhraw'
                call utmess('F', 'DVP_97', sk=saux08, si=codret)
            end if
            call as_mmhraw(idfimd, nomamd, nvtyge, atepai, nmatyp(ityp), &
                           zr(jr1eptu(ityp)), codret)
            if (codret .ne. 0) then
                saux08 = 'mmhraw'
                call utmess('F', 'DVP_97', sk=saux08, si=codret)
            end if
        end if
    end do cletype
!
    do ityp = 1, MT_NTYMAX
        if (nmatyp(ityp) .ne. 0) then
            call jedetr('&&IRMAES.CNX.'//nomtyp(ityp))
            if (exicoq .or. exigri) then
                call jedetr('&&IRMAES.EPAI.'//nomtyp(ityp))
            end if
            if (exipmf) then
                call jedetr('&&IRMAES.ORIP.'//nomtyp(ityp))
            end if
            if (exituy) then
                call jedetr('&&IRMAES.ORIT.'//nomtyp(ityp))
                call jedetr('&&IRMAES.RAYO.'//nomtyp(ityp))
            end if
        end if
    end do
!
end subroutine
