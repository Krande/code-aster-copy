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
subroutine rc36ac(noma, ncncin, chindi, chcara, nbma, &
                  listma, chresu)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/limend.h"
#include "asterfort/rc3601.h"
#include "asterfort/rc36sa.h"
#include "asterfort/rc36sn.h"
#include "asterfort/rc36sp.h"
#include "asterfort/rcma01.h"
#include "asterfort/rcmo01.h"
#include "asterfort/rcvale.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: nbma, listma(*)
    character(len=8) :: noma
    character(len=24) :: ncncin, chindi, chcara, chresu
!
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE_B3600
!
!     CALCUL DES AMPLITUDES DE CONTRAINTES
!     CALCUL DU FACTEUR D'USAGE
!
!     Pour chaque noeud de chaque maille:
!
!     pour une situation P, on a 2 états stabilisés
!     pour une situation Q, on a 2 états stabilisés
!
!     Soit 2 états stabilisés I et J appartenant respectivement aux
!     situations P et Q :
!
!     on calcule le SALT(I,J) = 0,5*(EC/E)*Ke*Sn(P,Q)*Sp(I,J)
!
!     avec Sn(P,Q) = Max( Sn(I,J) )
!          Sn(I,J) = Max( Max(Sn(I,J,ThP)), Max(Sn(I,J,ThQ)) )
!
!     avec Sp(I,J) = Max( Max(Sp(I,J,ThP)), Max(Sp(I,J,ThQ)) )
!
!
! Etape 1 : on calcule le SALT qui correspond aux combinaisons de tous
!           les états stabilisés appartenant aux situations d'un groupe
!           donné.
!
! Etape 2 : on calcule le SALT pour les situations non combinables
!
! Etape 3 : traitement des situations de passage
!           on calcule le SALT(I,J)
!              - avec I appartenant au premier groupe
!              - avec J appartenant au deuxieme groupe
!              - on lui associe le nombre d'occurrences de la
!                situation de passage
!
!
! IN  : NCNCIN : CONNECTIVITE INVERSE
! IN  : CHINDI : CHAM_ELEM DES INDICES DE CONTRAINTES
! IN  : CHCARA : CHAM_ELEM DES CARACTERISTIQUES ELEMENTAIRES
! IN  : NBMA   : NOMBRE DE MAILLES D'ANALYSE
! IN  : LISTMA : LISTE DES MAILLES D'ANALYSE
! OUT : CHRESU : CHAM_ELEM RESULTAT
!     ------------------------------------------------------------------
!
    integer(kind=8) :: ig, nbgr, nbsigr, jnsg, is1, ioc1, nocc, numgr, jcombi
    integer(kind=8) :: i1, nbth1, jth1, nbth2, nbcrs
    integer(kind=8) :: nbcin, nbcca, jcesl
    integer(kind=8) :: im, ima, nbpt, decrs, decin, decca, ipt, ino
    integer(kind=8) :: adrm, nbm, icmp, jconx1, jconx2, jfact, jpassa, npass
    integer(kind=8) :: ifm, niv, iocs, iad, ioc2, jcinl, jccal, nbp12, nbp23
    integer(kind=8) :: nbp13
    real(kind=8) :: ppi, ppj, snmax, samax, utot, saltij, ug, nadm(1), mpi(3)
    real(kind=8) :: mpj(3), sm, sn, sp, c(3), k(3), cara(3), matpi(14)
    real(kind=8) :: matpj(14), mse(3), snb, sab, smm, vale(2)
    aster_logical :: seisme, endur
    integer(kind=8) :: icodre(1)
    character(len=8) :: k8b, nommat, noeud, valk(7), kbid
    character(len=24) :: momepi, momepj, connex, matepi, matepj
    real(kind=8) :: typeke, spmeca, spther
    integer(kind=8), pointer :: situ_numero(:) => null()
    character(len=24), pointer :: materiau(:) => null()
    integer(kind=8), pointer :: situ_nb_occur(:) => null()
    real(kind=8), pointer :: situ_pres_a(:) => null()
    character(len=24), pointer :: situ_moment_b(:) => null()
    character(len=24), pointer :: situ_moment_a(:) => null()
    integer(kind=8), pointer :: ccad(:) => null()
    integer(kind=8), pointer :: cesd(:) => null()
    integer(kind=8), pointer :: cind(:) => null()
    real(kind=8), pointer :: situ_pres_b(:) => null()
    integer(kind=8), pointer :: situ_seisme(:) => null()
    real(kind=8), pointer :: ccav(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
    real(kind=8), pointer :: cinv(:) => null()
    character(len=8), pointer :: nom_materiau(:) => null()
    integer(kind=8), pointer :: situ_nume_group(:) => null()
! DEB ------------------------------------------------------------------
    call jemarq()
!
    call infniv(ifm, niv)
!
    connex = noma//'.CONNEX         '
    call jeveuo(connex, 'L', jconx1)
    call jeveuo(jexatr(connex, 'LONCUM'), 'L', jconx2)
!
    call jeveuo('&&RC3600.SITU_NUMERO', 'L', vi=situ_numero)
    call jelira('&&RC3600.SITU_NUME_GROUP', 'LONMAX', nbgr)
    call jeveuo('&&RC3600.SITU_NUME_GROUP', 'L', vi=situ_nume_group)
    call jeveuo('&&RC3600.SITU_SEISME', 'L', vi=situ_seisme)
!
    call jeveuo('&&RC3600.SITU_COMBINABLE', 'L', jcombi)
    call jeveuo('&&RC3600.SITU_PRES_A', 'L', vr=situ_pres_a)
    call jeveuo('&&RC3600.SITU_PRES_B', 'L', vr=situ_pres_b)
    call jeveuo('&&RC3600.SITU_MOMENT_A', 'L', vk24=situ_moment_a)
    call jeveuo('&&RC3600.SITU_MOMENT_B', 'L', vk24=situ_moment_b)
    call jeveuo('&&RC3600.SITU_NB_OCCUR', 'L', vi=situ_nb_occur)
    call jeveuo('&&RC3600.SITU_PASSAGE', 'L', jpassa)
!
    call jelira('&&RC32SI.PASSAGE_1_2', 'LONUTI', nbp12)
    call jelira('&&RC32SI.PASSAGE_2_3', 'LONUTI', nbp23)
    call jelira('&&RC32SI.PASSAGE_1_3', 'LONUTI', nbp13)
!
    call jeveuo('&&RC3600.MATERIAU', 'L', vk24=materiau)
    call jeveuo('&&RC3600.NOM_MATERIAU', 'L', vk8=nom_materiau)
!
! --- LE CHAM_ELEM RESULTAT
!
    call jeveuo(chresu(1:19)//'.CESD', 'L', vi=cesd)
    call jeveuo(chresu(1:19)//'.CESV', 'E', vr=cesv)
    call jeveuo(chresu(1:19)//'.CESL', 'E', jcesl)
    nbcrs = cesd(2)
!
! --- LE CHAMP INDICE DE CONTRAINTES
!
    call jeveuo(chindi(1:19)//'.CESV', 'L', vr=cinv)
    call jeveuo(chindi(1:19)//'.CESD', 'L', vi=cind)
    call jeveuo(chindi(1:19)//'.CESL', 'L', jcinl)
    nbcin = cind(2)
!
! --- LE CHAMP CARACTERISTIQUES
!
    call jeveuo(chcara(1:19)//'.CESV', 'L', vr=ccav)
    call jeveuo(chcara(1:19)//'.CESD', 'L', vi=ccad)
    call jeveuo(chcara(1:19)//'.CESL', 'L', jccal)
    nbcca = ccad(2)
!
    call wkvect('&&RC36AC_TRAVAIL', 'V V R', 4*50, jfact)
!
! --- IL FAUT CALCULER LE FACTEUR D'USAGE EN CHAQUE NOEUD DE CHAQUE
!     MAILLE
!
    do im = 1, nbma
!
        ima = listma(im)
        nommat = nom_materiau(ima)
!
        nbpt = cesd(5+4*(ima-1)+1)
        decrs = cesd(5+4*(ima-1)+4)
        decin = cind(5+4*(ima-1)+4)
        decca = ccad(5+4*(ima-1)+4)
!
        do ipt = 1, nbpt
!
! ------- LA CONNECTIVITE INVERSE
!
            ino = zi(jconx1-1+zi(jconx2+ima-1)+ipt-1)
            call jeveuo(jexnum(ncncin, ino), 'L', adrm)
            call jelira(jexnum(ncncin, ino), 'LONMAX', nbm)
            if (niv .ge. 2) then
                k8b = int_to_char8(ima)
                noeud = int_to_char8(ino)
                write (ifm, 1000) '===>> TRAITEMENT DU NOEUD ', noeud, &
                    ' APPARTENANT A LA MAILLE ', k8b
            end if
!
! ------- LES INDICES DE CONTRAINTES
!
            do icmp = 1, 3
                iad = decin+(ipt-1)*nbcin+icmp
                if (.not. zl(jcinl-1+iad)) then
                    valk(1) = int_to_char8(ino)
                    valk(2) = int_to_char8(ima)
                    if (icmp .eq. 1) then
                        valk(3) = 'C1'
                    else if (icmp .eq. 2) then
                        valk(3) = 'C2'
                    else
                        valk(3) = 'C3'
                    end if
                    call utmess('F', 'POSTRCCM_9', nk=3, valk=valk)
                end if
                c(icmp) = cinv(iad)
                iad = decin+(ipt-1)*nbcin+icmp+3
                if (.not. zl(jcinl-1+iad)) then
                    valk(1) = int_to_char8(ino)
                    valk(2) = int_to_char8(ima)
                    if (icmp .eq. 1) then
                        valk(3) = 'K1'
                    else if (icmp .eq. 2) then
                        valk(3) = 'K2'
                    else
                        valk(3) = 'K3'
                    end if
                    call utmess('F', 'POSTRCCM_9', nk=3, valk=valk)
                end if
                k(icmp) = cinv(iad)
            end do
!
! ------- LES CARATERISTIQUES : INERTIE, DIAMETRE, EPAISSEUR
!
            do icmp = 2, 4
                iad = decca+(ipt-1)*nbcca+icmp
                if (.not. zl(jccal-1+iad)) then
                    valk(1) = int_to_char8(ino)
                    valk(2) = int_to_char8(ima)
                    call utmess('F', 'POSTRCCM_8', nk=2, valk=valk)
                end if
                cara(icmp-1) = ccav(iad)
            end do
!
            sm = 0.d0
            snmax = 0.d0
            samax = 0.d0
            utot = 0.d0
!
! ----------------------------------------------------------------------
!                           E T A P E   1
! ----------------------------------------------------------------------
!
! ------- ON TRAITE LES SITUATIONS COMBINABLES DANS LEUR GROUPE
!         -----------------------------------------------------
!
            do ig = 1, nbgr
!
                numgr = situ_nume_group(ig)
                if (numgr .lt. 0) goto 100
                iocs = situ_seisme(ig)
!
                npass = 0
                seisme = .false.
!
                if (ig .eq. 1) then
                    if (nbp12 .ne. 0 .or. nbp13 .ne. 0) goto 100
                else if (ig .eq. 2) then
                    if (nbp12 .ne. 0 .or. nbp23 .ne. 0) goto 100
                else if (ig .eq. 3) then
                    if (nbp13 .ne. 0 .or. nbp23 .ne. 0) goto 100
                end if
!
! --------- PASSAGE 1 : PRISE EN COMPTE DU SEISME,
!                       CALCUL DU FACTEUR D'USAGE -> UTOT
!
                if (iocs .ne. 0) then
                    snb = 0.d0
                    sab = 0.d0
                    seisme = .true.
                    call rc3601(numgr, iocs, seisme, npass, ima, &
                                ipt, nbm, zi(adrm), c, k, &
                                cara, nommat, snb, sab, utot, &
                                sm, zr(jfact))
                    seisme = .false.
                end if
!
! --------- PASSAGE 2 : SANS LE SEISME
!                       CALCUL SU SN_MAX
!                       CALCUL SU SALT_MAX
!                       CALCUL DU FACTEUR D'USAGE -> UTOT
!
                call rc3601(numgr, iocs, seisme, npass, ima, &
                            ipt, nbm, zi(adrm), c, k, &
                            cara, nommat, snmax, samax, utot, &
                            sm, zr(jfact))
!
100             continue
            end do
!
! ----------------------------------------------------------------------
!                           E T A P E   2
! ----------------------------------------------------------------------
!
            mse(1) = 0.d0
            mse(2) = 0.d0
            mse(3) = 0.d0
!
! ------- ON TRAITE LES SITUATIONS NON COMBINABLES
!         ----------------------------------------
!
            do ig = 1, nbgr
!
                numgr = situ_nume_group(ig)
                if (numgr .lt. 0) goto 200
!
                call jelira(jexnum('&&RC3600.LES_GROUPES', numgr), 'LONMAX', nbsigr)
                call jeveuo(jexnum('&&RC3600.LES_GROUPES', numgr), 'L', jnsg)
!
                npass = 0
!
                do is1 = 1, nbsigr
                    ioc1 = zi(jnsg+is1-1)
                    if (zl(jcombi+ioc1-1)) goto 210
!
                    nocc = situ_nb_occur(1+2*ioc1-2)
!
                    ppi = situ_pres_a(ioc1)
                    momepi = situ_moment_a(ioc1)
                    call rcmo01(momepi, ima, ipt, mpi)
                    matepi = materiau(1+2*ioc1-1)
                    call rcma01(matepi, ima, ipt, nbm, zi(adrm), &
                                matpi)
!
                    ppj = situ_pres_b(ioc1)
                    momepj = situ_moment_b(ioc1)
                    call rcmo01(momepj, ima, ipt, mpj)
                    matepj = materiau(1+2*ioc1-2)
                    call rcma01(matepj, ima, ipt, nbm, zi(adrm), &
                                matpj)
!
                    call jelira(jexnum('&&RC3600.SITU_THERMIQUE', ioc1), 'LONUTI', nbth1)
                    if (nbth1 .ne. 0) then
                        call jeveuo(jexnum('&&RC3600.SITU_THERMIQUE', ioc1), 'L', jth1)
                    else
                        jth1 = 1
                    end if
!
                    nbth2 = 0
                    ioc2 = 0
!
! ----------- CALCUL DU SN
!
                    sn = 0.d0
                    call rc36sn(nbm, zi(adrm), ipt, c, cara, &
                                matpi, ppi, mpi, matpj, ppj, &
                                mpj, mse, nbth1, nbth2, ioc1, &
                                ioc2, sn)
                    snmax = max(snmax, sn)
!
! ----------- CALCUL DU SP
!
                    typeke = matpi(14)
                    sp = 0.d0
                    spmeca = 0.d0
                    spther = 0.d0
                    call rc36sp(nbm, zi(adrm), ipt, c, k, &
                                cara, matpi, ppi, mpi, matpj, &
                                ppj, mpj, mse, nbth1, nbth2, &
                                ioc1, ioc2, sp, typeke, spmeca, &
                                spther)
!
! ----------- CALCUL DU SALT
!
                    call rc36sa(nommat, matpi, matpj, sn, sp, &
                                typeke, spmeca, spther, saltij, smm)
!
                    if (saltij .gt. samax) then
                        samax = saltij
                        sm = smm
                    end if
!
! ----------- CALCUL DU FACTEUR D'USAGE
!
                    call limend(nommat, saltij, 'WOHLER', kbid, endur)
                    if (endur) then
                        ug = 0.d0
                    else
                        call rcvale(nommat, 'FATIGUE', 1, 'SIGM    ', [saltij], &
                                    1, 'WOHLER  ', nadm(1), icodre(1), 2)
                        if (nadm(1) .lt. 0) then
                            vale(1) = saltij
                            vale(2) = nadm(1)
                            call utmess('A', 'POSTRCCM_32', nr=2, valr=vale)
                        end if
                        ug = dble(nocc)/nadm(1)
                    end if
                    utot = utot+ug
!
210                 continue
                end do
!
200             continue
            end do
!
! ----------------------------------------------------------------------
!                           E T A P E   3
! ----------------------------------------------------------------------
!
! ------- ON TRAITE LES SITUATIONS DE PASSAGE
!         -----------------------------------
!
            do ig = 1, nbgr
!
                numgr = situ_nume_group(ig)
                if (numgr .ge. 0) goto 310
                numgr = -numgr
                iocs = situ_seisme(ig)
                if (iocs .eq. 0) then
                    seisme = .false.
                else
                    seisme = .true.
                end if
!
                call jelira(jexnum('&&RC3600.LES_GROUPES', numgr), 'LONMAX', nbsigr)
                call jeveuo(jexnum('&&RC3600.LES_GROUPES', numgr), 'L', jnsg)
                if (niv .ge. 2) then
                    write (ifm, 3004)
                    write (ifm, 3002) (situ_numero(1+zi(jnsg+i1-1)-1), i1=1, &
                                       nbsigr)
                end if
!
                npass = 7
!
                call rc3601(numgr, iocs, seisme, npass, ima, &
                            ipt, nbm, zi(adrm), c, k, &
                            cara, nommat, snmax, samax, utot, &
                            sm, zr(jfact))
!
310             continue
            end do
!
! ----------------------------------------------------------------------
!
! ------- ON STOCKE LES RESULTATS DE CALCUL
!         ---------------------------------
!
!         - LE SM
            icmp = 1
            iad = decrs+(ipt-1)*nbcrs+icmp
            cesv(iad) = sm
!         - LE SN
            icmp = 2
            iad = decrs+(ipt-1)*nbcrs+icmp
            cesv(iad) = snmax
!         - LE SN/3SM
            icmp = 3
            iad = decrs+(ipt-1)*nbcrs+icmp
            if (sm .eq. 0.d0) then
                cesv(iad) = 0.d0
            else
                cesv(iad) = snmax/(3*sm)
            end if
!         - LE SALT
            icmp = 4
            iad = decrs+(ipt-1)*nbcrs+icmp
            cesv(iad) = samax
!         - LE U_TOTAL
            icmp = 5
            iad = decrs+(ipt-1)*nbcrs+icmp
            cesv(iad) = utot
!
        end do
!
    end do
!
1000 format(a, a8, a, a8)
3002 format('=> LISTE DES NUMEROS DE SITUATION: ', 100(i4, 1x))
3004 format(/, '=> SITUATION DE PASSAGE')
!
    call jedema()
end subroutine
