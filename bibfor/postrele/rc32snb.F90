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
subroutine rc32snb(ze200, lieu, iocc1, iocc2, ns, &
                   sn, instsn, snet, sigmoypres, snther, &
                   sp3, spmeca3)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/rc32rt.h"
#include "asterfort/rc32s0.h"
#include "asterfort/rc32s0b.h"
#include "asterfort/rcZ2s0.h"
#include "asterfort/rctres.h"
#include "asterfort/utmess.h"
!
    character(len=4) :: lieu
    integer(kind=8) :: iocc1, iocc2, ns
    real(kind=8) :: sn, instsn(4), snet, sigmoypres, snther, sp3, spmeca3
    aster_logical :: ze200
!
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE ZE200 et B3200
!     CALCUL DU SN SANS SELECTION DES INSTANTS
!
!     ------------------------------------------------------------------
! IN  : ZE200  : EST ON EN ZE200 ou en B3200
! IN  : LIEU   : ='ORIG' : ORIGINE DU SEGMENT, ='EXTR' : EXTREMITE
! IN  : IOCC1  : NUMERO D'OCCURENCE DE LA PREMIERE SITUATION
! IN  : IOCC2  : NUMERO D'OCCURENCE DE LA DEUXIEME SITUATION
! IN  : NS     : 0 SANS SEISME OU 1 SI SEISME
! OUT : SN
! OUT : INSTSN : 2 INSTANTS DE SN ET 2 INSTANTS DE SN*
! OUT : SNET   : SN*
! OUT : SIGMOYPRES : CONTRAINTE MOYENNE DE PRESSION (pour le rochet)
! OUT : SNTHER : AMPLITUDE DE CONTRAINTES THERMIQUES(pour le rochet)
!
    integer(kind=8) :: jinfoi, nmecap, npresp, ntherp, jinfor, numcha, iret
    integer(kind=8) :: jchara, jcharb, k, j, jtranp, jsigu, i0, i1, e0(2)
    integer(kind=8) :: jtemp, i, l, nbinstp, nmecaq, npresq, ntherq, nbinstq
    integer(kind=8) :: jtranq, jvalin, jsnseis, inst1, inst2
    real(kind=8) :: presap, presbp, map(12), mbp(12), s2pp, sb(6), sbet(6)
    real(kind=8) :: smoypr1, instp(4)
    real(kind=8) :: tresca, sa(6), st(6), stet(6), smoypr(6), sc(6)
    real(kind=8) :: pij, tempa, tempb, A1p(12), B1p(12)
    real(kind=8) :: mom1(12), mom2(12), mij(12), sith(6), trescamax
    real(kind=8) :: trescaetmax, presaq, presbq, maq(12), mbq(12), sa1(6)
    real(kind=8) :: sa2(6), sa3(6), sa4(6), s2, st1(6), st2(6), st3(6)
    real(kind=8) :: A1q(12), B1q(12), st4(6), stet1(6), stet2(6), stet3(6)
    real(kind=8) :: stet4(6), momcstp(12), momcstq(12), k1, k2, k3, c1, c2
    real(kind=8) :: c3, snmec, snpres, sipr(6), simec(6), tresca1, tresca2
    real(kind=8) :: tresca3, tresca4
    character(len=8) :: knumec
    aster_logical :: tranp, tranq, unitaire
!
! DEB ------------------------------------------------------------------
    call jemarq()
!
    if (.not. ze200 .and. ns .ne. 0) call jeveuo('&&RC3200.SNSEISME.'//lieu, 'L', jsnseis)
!
    sn = 0.d0
    snet = 0.d0
    snther = 0.d0
    snmec = 0.d0
    snpres = 0.d0
    instsn(1) = -1.d0
    instsn(2) = -1.d0
    instsn(3) = -1.d0
    instsn(4) = -1.d0
    tranp = .false.
    tranq = .false.
    do k = 1, 12
        mom1(k) = 0.d0
        mom2(k) = 0.d0
        momcstp(k) = 0.d0
        momcstq(k) = 0.d0
    end do
    inst1 = 1
    inst2 = 1
!
! TROIS CONTRIBUTIONS POSSIBLES POUR SN
! SA  : UNITAIRE
!       SOUS SITUATION : CHAR_ETAT_A, CHAR_ETAT_B, PRES_A, PRES_B
! SB  : TRANSITOIRE
!       SOUS SITUATION : NUME_RESU_THER NUME_RESU_MECA, NUME_RESU_PRES
! SC  : INTERPOLATION DES MOMENTS
!       SOUS SITUATION : TEMP_A, TEMP_B, CHAR_ETAT_A, CHAR_ETAT_B, TABL_TEMP
!
!--------------------------------------------------------------------
!                  DANS LE CAS D'UNE SITUATION SEULE
!                          CALCUL DE SN(P,P)
!--------------------------------------------------------------------
    call jeveuo('&&RC3200.SITU_INFOI', 'L', jinfoi)
    call jeveuo('&&RC3200.SITU_INFOR', 'L', jinfor)
!
!-- on regarde si la pression est sous forme unitaire ou transitoire
!-- on regarde si la méca est sous forme unitaire ou transitoire
    nmecap = zi(jinfoi+27*(iocc1-1)+23)
    npresp = zi(jinfoi+27*(iocc1-1)+22)
    ntherp = zi(jinfoi+27*(iocc1-1)+26)
!
!---- Mot clé pour accélerer les calculs si aucun chargement en unitaire
    unitaire = .false.
    if (npresp .eq. 1 .or. nmecap .eq. 1) unitaire = .true.
!
    presap = zr(jinfor+4*(iocc1-1))
    presbp = zr(jinfor+4*(iocc1-1)+1)
!
    do k = 1, 12
        map(k) = 0.d0
        mbp(k) = 0.d0
    end do
!
    if (nmecap .eq. 1 .or. nmecap .eq. 3) then
!------ Chargement état A
        numcha = zi(jinfoi+27*(iocc1-1)+24)
        knumec = 'C       '
        call codent(numcha, 'D0', knumec(2:8))
        call jeexin(jexnom('&&RC3200.VALE_CHAR', knumec), iret)
        if (iret .eq. 0) call utmess('F', 'POSTRCCM_51')
        call jeveuo(jexnom('&&RC3200.VALE_CHAR', knumec), 'L', jchara)
!------ Chargement état B
        numcha = zi(jinfoi+27*(iocc1-1)+25)
        knumec = 'C       '
        call codent(numcha, 'D0', knumec(2:8))
        call jeexin(jexnom('&&RC3200.VALE_CHAR', knumec), iret)
        if (iret .eq. 0) call utmess('F', 'POSTRCCM_51')
        call jeveuo(jexnom('&&RC3200.VALE_CHAR', knumec), 'L', jcharb)
!
        do k = 1, 12
            map(k) = zr(jchara-1+k)
            mbp(k) = zr(jcharb-1+k)
        end do
!
    end if
!-- On récupère les contraintes unitaires
    if (.not. ze200) then
        if (nmecap .eq. 1 .or. nmecap .eq. 3 .or. npresp .eq. 1) then
            call jeveuo('&&RC3200.MECA_UNIT .'//lieu, 'L', jsigu)
        end if
    end if
!
!-- SA partie unitaire (moments et/ou pression)
    do j = 1, 6
        sa(j) = 0.d0
    end do
    if (.not. ze200 .and. nmecap .eq. 1) then
        do j = 1, 6
            do k = 1, 12
                sa(j) = sa(j)+(map(k)-mbp(k))*zr(jsigu-1+78+6*(k-1)+j)
            end do
        end do
    end if
    if (.not. ze200 .and. npresp .eq. 1) then
        do j = 1, 6
            sa(j) = sa(j)+(presap-presbp)*zr(jsigu-1+78+72+j)
        end do
    end if
!
!-- SB partie transitoire
    do j = 1, 6
        sb(j) = 0.d0
        sbet(j) = 0.d0
        st(j) = 0.d0
        stet(j) = 0.d0
        smoypr(j) = 0.d0
    end do
!
    if (nmecap .eq. 2 .or. npresp .eq. 2 .or. ntherp .eq. 1) then
        call jeveuo(jexnum('&&RC3200.TRANSIT.'//lieu, iocc1), 'L', jtranp)
        tranp = .true.
        nbinstp = int(zr(jtranp))
    else
        instp(1) = -1.0
        instp(2) = -1.0
        instp(3) = -1.0
        instp(4) = -1.0
    end if
!
!-- SC partie unitaire avec interpolation des moments
    do j = 1, 6
        sc(j) = 0.d0
    end do
!
    if (nmecap .eq. 3) then
        tempa = zr(jinfor+4*(iocc1-1)+2)
        tempb = zr(jinfor+4*(iocc1-1)+3)
        if (tempa .lt. tempb) then
            do k = 1, 12
                A1p(k) = (map(k)-mbp(k))/(tempa-tempb)
                B1p(k) = (mbp(k)*tempa-map(k)*tempb)/(tempa-tempb)
            end do
        else
            do k = 1, 12
                A1p(k) = (mbp(k)-map(k))/(tempb-tempa)
                B1p(k) = (map(k)*tempb-mbp(k)*tempa)/(tempb-tempa)
            end do
        end if
!
! ----  si la situation ne possède aucun transitoire
! ----  (ni thermique, ni de pression ni mécanique)
        if (nmecap .ne. 2 .and. npresp .ne. 2 .and. ntherp .ne. 1) then
            call jeveuo(jexnum('&&RC3200.TEMPCST', iocc1), 'L', jtemp)
            do k = 1, 12
                if (lieu .eq. 'ORIG') then
                    momcstp(k) = A1p(k)*zr(jtemp)+B1p(k)
                else
                    momcstp(k) = A1p(k)*zr(jtemp+1)+B1p(k)
                end if
            end do
            do j = 1, 6
                do k = 1, 12
                    sc(j) = sc(j)+momcstp(k)*zr(jsigu-1+78+6*(k-1)+j)
                end do
            end do
        end if
    end if
!--------------------------------------------------------------------
!                  DANS LE CAS D'UNE SITUATION SEULE
!--------------------------------------------------------------------
!-- CALCUL DU SN selon le chapitre du RCC-M étudié
!
    do i0 = 1, 2
        i1 = 2*(i0-2)+1
        e0(i0) = i1
    end do
!
    trescamax = -1.d0
    trescaetmax = -1.d0
    smoypr1 = 0.d0
!
    if (iocc2 .eq. 0 .or. iocc1 .eq. iocc2) then
! ----- si on est en ZE200
        if (ze200) then
            call rcZ2s0('SN', map, mbp, presap, presbp, &
                        ns, s2pp)
            if (tranp) then
                do i = 1, nbinstp
                    do l = 1, nbinstp
                        do j = 1, 6
                            sb(j) = zr(jtranp+50*(i-1)+1+12+j)-zr(jtranp+50*(l-1)+1+12+j)
                            sbet(j) = ( &
                                      zr( &
                                      jtranp+50*(i-1)+1+12+j)-zr(jtranp+50*(i-1)+1+42+j))-(zr(j&
                                      &tranp+50*(l-1)+1+12+j)-zr(jtranp+50*(l-1)+1+42+j &
                                      ) &
                                      )
                        end do
                        call rctres(sb, tresca)
                        if (tresca .gt. trescamax) then
                            trescamax = tresca
                            instp(1) = zr(jtranp+50*(i-1)+1)
                            instp(2) = zr(jtranp+50*(l-1)+1)
                            inst1 = i
                            inst2 = l
                        end if
                        call rctres(sbet, tresca)
                        if (tresca .gt. trescaetmax) then
                            trescaetmax = tresca
                            instp(3) = zr(jtranp+50*(i-1)+1)
                            instp(4) = zr(jtranp+50*(l-1)+1)
                        end if
                    end do
                end do
                sn = s2pp+trescamax
                snet = s2pp+trescaetmax
                do j = 1, 6
                    sith(j) = zr(jtranp+50*(inst1-1)+1+24+j)-zr(jtranp+50*(inst2-1)+1+24+j)
                    sipr(j) = zr(jtranp+50*(inst1-1)+1+30+j)-zr(jtranp+50*(inst2-1)+1+30+j)
                end do
                call rctres(sith, snther)
                call rctres(sipr, snpres)
            else
                sn = s2pp
                snet = s2pp
            end if
        else
! ----- si on est en B3200
            if (tranp) then
                do i = 1, nbinstp
                    do l = 1, nbinstp
!
                        do i0 = 1, 2
                            do j = 1, 6
                                sc(j) = 0.d0
                                sb(j) = zr(jtranp+50*(i-1)+1+12+j)-zr(jtranp+50*(l-1)+1+12+j)
                                sbet(j) = ( &
                                          zr( &
                                          jtranp+50*(i-1)+1+12+j)-zr(jtranp+50*(i-1)+1+42+j))-(&
                                          &zr(jtranp+50*(l-1)+1+12+j)-zr(jtranp+50*(l-1)+1+42+j &
                                          ) &
                                          )
!
                                if (nmecap .eq. 3) then
                                    do k = 1, 12
                             mij(k) = A1p(k)*(zr(jtranp+50*(i-1)+1+48+1)-zr(jtranp+50*(l-1)+1+48+1))
                                        sc(j) = sc(j)+mij(k)*zr(jsigu-1+78+6*(k-1)+j)
                                    end do
                                end if
!
                                st(j) = sb(j)+sc(j)+e0(i0)*sa(j)
                                stet(j) = sbet(j)+sc(j)+e0(i0)*sa(j)
                                smoypr(j) = zr(jtranp+50*(i-1)+1+18+j)
                            end do
                            call rctres(smoypr, tresca)
                            smoypr1 = max(smoypr1, tresca)
                            if (ns .eq. 0) then
                                call rctres(st, tresca)
                                if (tresca .gt. trescamax) then
                                    trescamax = tresca
                                    inst1 = i
                                    inst2 = l
                                    instp(1) = zr(jtranp+50*(i-1)+1)
                                    instp(2) = zr(jtranp+50*(l-1)+1)
                                end if
                                call rctres(stet, tresca)
                                if (tresca .gt. trescaetmax) then
                                    trescaetmax = tresca
                                    instp(3) = zr(jtranp+50*(i-1)+1)
                                    instp(4) = zr(jtranp+50*(l-1)+1)
                                end if
                            else
                                call rc32s0b(zr(jsnseis), st, tresca)
                                if (tresca .gt. trescamax) then
                                    trescamax = tresca
                                    instp(1) = zr(jtranp+50*(i-1)+1)
                                    instp(2) = zr(jtranp+50*(l-1)+1)
                                    inst1 = i
                                    inst2 = l
                                end if
                                call rc32s0b(zr(jsnseis), stet, tresca)
                                if (tresca .gt. trescaetmax) then
                                    trescaetmax = tresca
                                    instp(3) = zr(jtranp+50*(i-1)+1)
                                    instp(4) = zr(jtranp+50*(l-1)+1)
                                end if
                            end if
                            if (.not. unitaire) exit
                        end do
                    end do
                end do
!
                do j = 1, 6
                    sith(j) = zr(jtranp+50*(inst1-1)+1+24+j)-zr(jtranp+50*(inst2-1)+1+24+j)
                    sipr(j) = zr(jtranp+50*(inst1-1)+1+30+j)-zr(jtranp+50*(inst2-1)+1+30+j)
                    simec(j) = zr(jtranp+50*(inst1-1)+1+36+j)-zr(jtranp+50*(inst2-1)+1+36+j)
                end do
                call rctres(sith, snther)
                call rctres(sipr, snpres)
                if (ns .eq. 0) then
                    call rctres(simec, snmec)
                else
                    call rc32s0b(zr(jsnseis), simec, snmec)
                end if
!
            else
                instp(1) = -1.d0
                instp(2) = -1.d0
                instp(3) = -1.d0
                instp(4) = -1.d0
                do i0 = 1, 2
                    do j = 1, 6
                        st(j) = sc(j)+e0(i0)*sa(j)
                        stet(j) = sc(j)+e0(i0)*sa(j)
                    end do
                    if (ns .eq. 0) then
                        call rctres(st, tresca)
                        trescamax = max(tresca, trescamax)
                        call rctres(stet, tresca)
                        trescaetmax = max(tresca, trescaetmax)
                    else
                        call rc32s0b(zr(jsnseis), st, tresca)
                        trescamax = max(tresca, trescamax)
                        call rc32s0b(zr(jsnseis), stet, tresca)
                        trescaetmax = max(tresca, trescaetmax)
                    end if
                end do
            end if
            sn = trescamax
            snet = trescaetmax
        end if
!
        instsn(1) = instp(1)
        instsn(2) = instp(2)
        instsn(3) = instp(3)
        instsn(4) = instp(4)
!
!-- Pour le calcul du rochet thermique
        sigmoypres = r8vide()
        do j = 1, 6
            smoypr(j) = 0.d0
        end do
        if (npresp .eq. 1) then
            if (ze200) then
                call rc32rt(presap, presbp, sigmoypres)
            else
                do j = 1, 6
                    pij = max(abs(presap), abs(presbp))
                    smoypr(j) = pij*zr(jsigu-1+156+72+j)
                end do
                call rctres(smoypr, sigmoypres)
            end if
        else if (npresp .eq. 2) then
            sigmoypres = smoypr1
        end if
!
    else
!--------------------------------------------------------------------
!                  DANS LE CAS D'UNE COMBINAISON DE SITUATIONS
!--------------------------------------------------------------------
!-- on regarde si la pression est sous forme unitaire ou transitoire
!-- on regarde si la méca est sous forme unitaire ou transitoire
        nmecaq = zi(jinfoi+27*(iocc2-1)+23)
        npresq = zi(jinfoi+27*(iocc2-1)+22)
        ntherq = zi(jinfoi+27*(iocc2-1)+26)
!
        presaq = zr(jinfor+4*(iocc2-1))
        presbq = zr(jinfor+4*(iocc2-1)+1)
!
        do k = 1, 12
            maq(k) = 0.d0
            mbq(k) = 0.d0
        end do
!
!-- On récupère les contraintes unitaires
        if (.not. ze200) then
            if (nmecaq .eq. 1 .or. nmecaq .eq. 3 .or. npresq .eq. 1) then
                call jeveuo('&&RC3200.MECA_UNIT .'//lieu, 'L', jsigu)
            end if
        end if
!
        if (nmecaq .eq. 1 .or. nmecaq .eq. 3) then
!------ Chargement état A
            numcha = zi(jinfoi+27*(iocc2-1)+24)
            knumec = 'C       '
            call codent(numcha, 'D0', knumec(2:8))
            call jeexin(jexnom('&&RC3200.VALE_CHAR', knumec), iret)
            if (iret .eq. 0) call utmess('F', 'POSTRCCM_51')
            call jeveuo(jexnom('&&RC3200.VALE_CHAR', knumec), 'L', jchara)
!------ Chargement état B
            numcha = zi(jinfoi+27*(iocc2-1)+25)
            knumec = 'C       '
            call codent(numcha, 'D0', knumec(2:8))
            call jeexin(jexnom('&&RC3200.VALE_CHAR', knumec), iret)
            if (iret .eq. 0) call utmess('F', 'POSTRCCM_51')
            call jeveuo(jexnom('&&RC3200.VALE_CHAR', knumec), 'L', jcharb)
            do k = 1, 12
                maq(k) = zr(jchara-1+k)
                mbq(k) = zr(jcharb-1+k)
            end do
        end if
!
!-- SA partie unitaire (moments et/ou pression)
        do j = 1, 6
            sa1(j) = 0.d0
            sa2(j) = 0.d0
            sa3(j) = 0.d0
            sa4(j) = 0.d0
        end do
        if (.not. ze200) then
            if (nmecap .eq. 1) then
                do j = 1, 6
                    do k = 1, 12
                        sa1(j) = sa1(j)+map(k)*zr(jsigu-1+78+6*(k-1)+j)
                        sa2(j) = sa2(j)+map(k)*zr(jsigu-1+78+6*(k-1)+j)
                        sa3(j) = sa3(j)+mbp(k)*zr(jsigu-1+78+6*(k-1)+j)
                        sa4(j) = sa4(j)+mbp(k)*zr(jsigu-1+78+6*(k-1)+j)
                    end do
                end do
            end if
!
            if (nmecaq .eq. 1) then
                do j = 1, 6
                    do k = 1, 12
                        sa1(j) = sa1(j)-mbq(k)*zr(jsigu-1+78+6*(k-1)+j)
                        sa2(j) = sa2(j)-maq(k)*zr(jsigu-1+78+6*(k-1)+j)
                        sa3(j) = sa3(j)-mbq(k)*zr(jsigu-1+78+6*(k-1)+j)
                        sa4(j) = sa4(j)-maq(k)*zr(jsigu-1+78+6*(k-1)+j)
                    end do
                end do
            end if
!
            if (npresp .eq. 1 .or. npresq .eq. 1) then
                do j = 1, 6
                    sa1(j) = sa1(j)+(presap-presbq)*zr(jsigu-1+78+72+j)
                    sa2(j) = sa2(j)+(presap-presaq)*zr(jsigu-1+78+72+j)
                    sa3(j) = sa3(j)+(presbp-presbq)*zr(jsigu-1+78+72+j)
                    sa4(j) = sa4(j)+(presbp-presaq)*zr(jsigu-1+78+72+j)
                end do
            end if
        end if
!
!-- SB partie transitoire
        do j = 1, 6
            sb(j) = 0.d0
            sbet(j) = 0.d0
            st1(j) = 0.d0
            st2(j) = 0.d0
            st3(j) = 0.d0
            st4(j) = 0.d0
            stet1(j) = 0.d0
            stet2(j) = 0.d0
            stet3(j) = 0.d0
            stet4(j) = 0.d0
        end do
!
        if (nmecaq .eq. 2 .or. npresq .eq. 2 .or. ntherq .eq. 1) then
            call jeveuo(jexnum('&&RC3200.TRANSIT.'//lieu, iocc2), 'L', jtranq)
            tranq = .true.
            nbinstq = int(zr(jtranq))
        end if
!
!-- SC partie unitaire avec interpolation des moments
        do j = 1, 6
            sc(j) = 0.d0
        end do
!
        if (nmecaq .eq. 3) then
            tempa = zr(jinfor+4*(iocc2-1)+2)
            tempb = zr(jinfor+4*(iocc2-1)+3)
            if (tempa .lt. tempb) then
                do k = 1, 12
                    A1q(k) = (maq(k)-mbq(k))/(tempa-tempb)
                    B1q(k) = (mbq(k)*tempa-maq(k)*tempb)/(tempa-tempb)
                end do
            else
                do k = 1, 12
                    A1q(k) = (mbq(k)-maq(k))/(tempb-tempa)
                    B1q(k) = (maq(k)*tempb-mbq(k)*tempa)/(tempb-tempa)
                end do
            end if
!
! ----  si la situation ne possède aucun transitoire
! ----  (ni thermique, ni de pression ni mécanique)
            if (.not. tranq) then
                call jeveuo(jexnum('&&RC3200.TEMPCST', iocc2), 'L', jtemp)
                do k = 1, 12
                    if (lieu .eq. 'ORIG') then
                        momcstq(k) = A1q(k)*zr(jtemp)+B1q(k)
                    else
                        momcstq(k) = A1q(k)*zr(jtemp+1)+B1q(k)
                    end if
                end do
            end if
        end if
!
! --------------------------------------------------------------
!                          CALCUL DE SN(P,Q)
! --------------------------------------------------------------
!
        trescamax = -1.d0
        trescaetmax = -1.d0
        tresca1 = 0.d0
        tresca2 = 0.d0
        tresca3 = 0.d0
        tresca4 = 0.d0
        s2 = 0.d0
        unitaire = .false.
!
!---- Mot clé pour accélerer les calculs si aucun chargement en unitaire
        if (npresp .eq. 1 .or. nmecap .eq. 1) unitaire = .true.
        if (npresq .eq. 1 .or. nmecaq .eq. 1) unitaire = .true.
!
        if (ze200) then
! ----- si on est en ZE200
            call rcZ2s0('SN', map, mbq, presap, presbq, &
                        ns, s2pp)
            s2 = max(s2, s2pp)
            call rcZ2s0('SN', map, maq, presap, presaq, &
                        ns, s2pp)
            s2 = max(s2, s2pp)
            call rcZ2s0('SN', mbp, maq, presbp, presaq, &
                        ns, s2pp)
            s2 = max(s2, s2pp)
            call rcZ2s0('SN', mbp, mbq, presbp, presbq, &
                        ns, s2pp)
            s2 = max(s2, s2pp)
!
            if (tranp) then
                if (tranq) then
                    do i = 1, nbinstp
                        do l = 1, nbinstq
                            do j = 1, 6
                                sb(j) = zr(jtranp+50*(i-1)+1+12+j)-zr(jtranq+50*(l-1)+1+12+j)
                                sbet(j) = ( &
                                          zr( &
                                          jtranp+50*(i-1)+1+12+j)-zr(jtranp+50*(i-1)+1+42+j))-(&
                                          &zr(jtranq+50*(l-1)+1+12+j)-zr(jtranq+50*(l-1)+1+42+j &
                                          ) &
                                          )
                            end do
                            call rctres(sb, tresca)
                            if (tresca .gt. trescamax) then
                                trescamax = tresca
                                inst1 = i
                                inst2 = l
                            end if
                            call rctres(sbet, tresca)
                            if (tresca .gt. trescaetmax) then
                                trescaetmax = tresca
                                instp(3) = zr(jtranp+50*(i-1)+1)
                                instp(4) = zr(jtranq+50*(l-1)+1)
                            end if
                        end do
                    end do
!
                    instp(1) = zr(jtranp+50*(inst1-1)+1)
                    instp(2) = zr(jtranq+50*(inst2-1)+1)
                    do j = 1, 6
                        sith(j) = zr(jtranp+50*(inst1-1)+1+24+j)-zr(jtranq+50*(inst2-1)+1+24+j)
                        sipr(j) = zr(jtranp+50*(inst1-1)+1+30+j)-zr(jtranq+50*(inst2-1)+1+30+j)
                    end do
                    call rctres(sith, snther)
                    call rctres(sipr, snpres)
!
                else
                    do i = 1, nbinstp
                        do j = 1, 6
                            sb(j) = zr(jtranp+50*(i-1)+1+12+j)
                            sbet(j) = (zr(jtranp+50*(i-1)+1+12+j)-zr(jtranp+50*(i-1)+1+42+j))
                        end do
                        call rctres(sb, tresca)
                        if (tresca .gt. trescamax) then
                            trescamax = tresca
                            inst1 = i
                        end if
                        call rctres(sbet, tresca)
                        if (tresca .gt. trescaetmax) then
                            trescaetmax = tresca
                            instp(3) = zr(jtranp+50*(i-1)+1)
                            instp(4) = -1.d0
                        end if
                    end do
!
                    instp(1) = zr(jtranp+50*(inst1-1)+1)
                    instp(2) = -1.d0
                    do j = 1, 6
                        sith(j) = zr(jtranp+50*(inst1-1)+1+24+j)
                        sipr(j) = zr(jtranp+50*(inst1-1)+1+30+j)
                    end do
                    call rctres(sith, snther)
                    call rctres(sipr, snpres)
!
                end if
            else
                if (tranq) then
                    do l = 1, nbinstq
                        do j = 1, 6
                            sb(j) = -zr(jtranq+50*(l-1)+1+12+j)
                            sbet(j) = -(zr(jtranq+50*(l-1)+1+12+j)-zr(jtranq+50*(l-1)+1+42+j))
                        end do
                        call rctres(sb, tresca)
                        if (tresca .gt. trescamax) then
                            trescamax = tresca
                            inst2 = l
                        end if
                        call rctres(sbet, tresca)
                        if (tresca .gt. trescaetmax) then
                            trescaetmax = tresca
                            instp(3) = -1.d0
                            instp(4) = zr(jtranq+50*(l-1)+1)
                        end if
                    end do
!
                    instp(1) = -1.d0
                    instp(2) = zr(jtranq+50*(inst2-1)+1)
                    do j = 1, 6
                        sith(j) = -zr(jtranq+50*(inst2-1)+1+24+j)
                        sipr(j) = -zr(jtranq+50*(inst2-1)+1+30+j)
                    end do
                    call rctres(sith, snther)
                    call rctres(sipr, snpres)
!
                else
                    instp(1) = -1.d0
                    instp(2) = -1.d0
                    instp(3) = -1.d0
                    instp(4) = -1.d0
                end if
!
            end if
!
            sn = s2+trescamax
            snet = s2+trescaetmax
!
        else
! ----- si on est en B3200
!
            if (tranp) then
! ----- la première situation a un transitoire
                if (tranq) then
! ------- la deuxième situation a un transitoire
                    do i = 1, nbinstp
                        do l = 1, nbinstq
                            do i0 = 1, 2
                                do j = 1, 6
                                    sc(j) = 0.d0
                                    sb(j) = zr(jtranp+50*(i-1)+1+12+j)-zr(jtranq+50*(l-1)+1+12+j &
                                                                          )
                                    sbet(j) = ( &
                                              zr( &
                                              jtranp+50*(i-1)+1+12+j)-zr(jtranp+50*(i-1)+1+42+j))&
                                              & -(zr(jtranq+50*(l-1)+1+12+j)-zr(jtranq+50*(l-1)+&
                                              &1+42+j &
                                              ) &
                                              )
                                    do k = 1, 12
                                        if (nmecap .eq. 3) mom1(k) = A1p(k)*zr(jtranp+50*(i-1)+1+&
                                                                     &48+1)+B1p(k)
                                        if (nmecaq .eq. 3) mom2(k) = A1q(k)*zr(jtranq+50*(l-1)+1+&
                                                                     &48+1)+B1q(k)
                                        mij(k) = mom1(k)-mom2(k)
                                        if (nmecap .eq. 3 .or. nmecaq .eq. 3) then
                                            sc(j) = sc(j)+mij(k)*zr(jsigu-1+78+6*(k-1)+j)
                                        end if
                                    end do
                                    st1(j) = sb(j)+sc(j)+e0(i0)*sa1(j)
                                    stet1(j) = sbet(j)+sc(j)+e0(i0)*sa1(j)
                                    if (unitaire) then
                                        st2(j) = sb(j)+sc(j)+e0(i0)*sa2(j)
                                        st3(j) = sb(j)+sc(j)+e0(i0)*sa3(j)
                                        st4(j) = sb(j)+sc(j)+e0(i0)*sa4(j)
                                        stet2(j) = sbet(j)+sc(j)+e0(i0)*sa2(j)
                                        stet3(j) = sbet(j)+sc(j)+e0(i0)*sa3(j)
                                        stet4(j) = sbet(j)+sc(j)+e0(i0)*sa4(j)
                                    end if
                                end do
!
!------------------ Calcul du SN sans séisme
                                if (ns .eq. 0) then
                                    call rctres(st1, tresca1)
                                    if (unitaire) then
                                        call rctres(st2, tresca2)
                                        call rctres(st3, tresca3)
                                        call rctres(st4, tresca4)
                                    end if
                                    tresca = max(tresca1, tresca2, tresca3, tresca4)
!
                                    if (tresca .gt. trescamax) then
                                        trescamax = tresca
                                        inst1 = i
                                        inst2 = l
                                    end if
!------------------ Calcul du SN avec séisme
                                else
                                    call rc32s0b(zr(jsnseis), st1, tresca1)
                                    if (unitaire) then
                                        call rc32s0b(zr(jsnseis), st2, tresca2)
                                        call rc32s0b(zr(jsnseis), st3, tresca3)
                                        call rc32s0b(zr(jsnseis), st4, tresca4)
                                    end if
                                    tresca = max(tresca1, tresca2, tresca3, tresca4)
!
                                    if (tresca .gt. trescamax) then
                                        trescamax = tresca
                                        inst1 = i
                                        inst2 = l
                                    end if
                                end if
!------------------ Calcul du SN* sans séisme
                                if (ns .eq. 0) then
                                    call rctres(stet1, tresca1)
                                    if (unitaire) then
                                        call rctres(stet2, tresca2)
                                        call rctres(stet3, tresca3)
                                        call rctres(stet4, tresca4)
                                    end if
                                    tresca = max(tresca1, tresca2, tresca3, tresca4)
                                    if (tresca .gt. trescaetmax) then
                                        trescaetmax = tresca
                                        instp(3) = zr(jtranp+50*(i-1)+1)
                                        instp(4) = zr(jtranq+50*(l-1)+1)
                                    end if
                                end if
!------------------ Calcul du SN* avec séisme
                                if (ns .ne. 0) then
                                    call rc32s0b(zr(jsnseis), stet1, tresca1)
                                    if (unitaire) then
                                        call rc32s0b(zr(jsnseis), stet2, tresca2)
                                        call rc32s0b(zr(jsnseis), stet3, tresca3)
                                        call rc32s0b(zr(jsnseis), stet4, tresca4)
                                    end if
                                    tresca = max(tresca1, tresca2, tresca3, tresca4)
                                    if (tresca .gt. trescaetmax) then
                                        trescaetmax = tresca
                                        instp(3) = zr(jtranp+50*(i-1)+1)
                                        instp(4) = zr(jtranq+50*(l-1)+1)
                                    end if
                                end if
                                if (.not. unitaire) exit
!
                            end do
                        end do
                    end do
!
!---------- Calcul de snther, snpres et snmec
                    do j = 1, 6
                        sith(j) = zr(jtranp+50*(inst1-1)+1+24+j)-zr(jtranq+50*(inst2-1)+1+24+j)
                        sipr(j) = zr(jtranp+50*(inst1-1)+1+30+j)-zr(jtranq+50*(inst2-1)+1+30+j)
                        simec(j) = zr(jtranp+50*(inst1-1)+1+36+j)-zr(jtranq+50*(inst2-1)+1+36+j)
                    end do
                    call rctres(sith, snther)
                    call rctres(sipr, snpres)
                    if (ns .eq. 0) then
                        call rctres(simec, snmec)
                    else
                        call rc32s0b(zr(jsnseis), simec, snmec)
                    end if
                    instp(1) = zr(jtranp+50*(inst1-1)+1)
                    instp(2) = zr(jtranq+50*(inst2-1)+1)
!
                else
! ------- la deuxième situation n'a pas de transitoire
                    do i = 1, nbinstp
                        do i0 = 1, 2
                            do j = 1, 6
                                sc(j) = 0.d0
                                sb(j) = zr(jtranp+50*(i-1)+1+12+j)-0.d0
                                sbet(j) = (zr(jtranp+50*(i-1)+1+12+j)-zr(jtranp+50*(i-1)+1+42+j))
                                do k = 1, 12
                                    if (nmecap .eq. 3) mom1(k) = A1p(k)*zr(jtranp+50*(i-1)+1+48+1&
                                                                 &)+B1p(k)
                                    mij(k) = mom1(k)-momcstq(k)
                                    if (nmecap .eq. 3 .or. nmecaq .eq. 3) then
                                        sc(j) = sc(j)+mij(k)*zr(jsigu-1+78+6*(k-1)+j)
                                    end if
                                end do
                                st1(j) = sb(j)+sc(j)+e0(i0)*sa1(j)
                                stet1(j) = sbet(j)+sc(j)+e0(i0)*sa1(j)
                                if (unitaire) then
                                    st2(j) = sb(j)+sc(j)+e0(i0)*sa2(j)
                                    st3(j) = sb(j)+sc(j)+e0(i0)*sa3(j)
                                    st4(j) = sb(j)+sc(j)+e0(i0)*sa4(j)
                                    stet2(j) = sbet(j)+sc(j)+e0(i0)*sa2(j)
                                    stet3(j) = sbet(j)+sc(j)+e0(i0)*sa3(j)
                                    stet4(j) = sbet(j)+sc(j)+e0(i0)*sa4(j)
                                end if
                            end do
!
!------------------ Calcul du SN sans séisme
                            if (ns .eq. 0) then
                                call rctres(st1, tresca1)
                                if (unitaire) then
                                    call rctres(st2, tresca2)
                                    call rctres(st3, tresca3)
                                    call rctres(st4, tresca4)
                                end if
                                tresca = max(tresca1, tresca2, tresca3, tresca4)
!
                                if (tresca .gt. trescamax) then
                                    trescamax = tresca
                                    inst1 = i
                                end if
                            end if
!------------------ Calcul du SN avec séisme
                            if (ns .ne. 0) then
                                call rc32s0b(zr(jsnseis), st1, tresca1)
                                if (unitaire) then
                                    call rc32s0b(zr(jsnseis), st2, tresca2)
                                    call rc32s0b(zr(jsnseis), st3, tresca3)
                                    call rc32s0b(zr(jsnseis), st4, tresca4)
                                end if
                                tresca = max(tresca1, tresca2, tresca3, tresca4)
!
                                if (tresca .gt. trescamax) then
                                    trescamax = tresca
                                    inst1 = i
                                end if
                            end if
!------------------ Calcul du SN* sans séisme
                            if (ns .eq. 0) then
                                call rctres(stet1, tresca1)
                                if (unitaire) then
                                    call rctres(stet2, tresca2)
                                    call rctres(stet3, tresca3)
                                    call rctres(stet4, tresca4)
                                end if
                                tresca = max(tresca1, tresca2, tresca3, tresca4)
                                if (tresca .gt. trescaetmax) then
                                    trescaetmax = tresca
                                    instp(3) = zr(jtranp+50*(i-1)+1)
                                    instp(4) = -1.d0
                                end if
                            end if
!------------------ Calcul du SN* avec séisme
                            if (ns .ne. 0) then
                                call rc32s0b(zr(jsnseis), stet1, tresca1)
                                if (unitaire) then
                                    call rc32s0b(zr(jsnseis), stet2, tresca2)
                                    call rc32s0b(zr(jsnseis), stet3, tresca3)
                                    call rc32s0b(zr(jsnseis), stet4, tresca4)
                                end if
                                tresca = max(tresca1, tresca2, tresca3, tresca4)
                                if (tresca .gt. trescaetmax) then
                                    trescaetmax = tresca
                                    instp(3) = zr(jtranp+50*(i-1)+1)
                                    instp(4) = -1.d0
                                end if
                            end if
                            if (.not. unitaire) exit
!
                        end do
                    end do
!
!---------- Calcul de snther, snpres et snmec
                    do j = 1, 6
                        sith(j) = zr(jtranp+50*(inst1-1)+1+24+j)
                        sipr(j) = zr(jtranp+50*(inst1-1)+1+30+j)
                        simec(j) = zr(jtranp+50*(inst1-1)+1+36+j)
                    end do
                    call rctres(sith, snther)
                    call rctres(sipr, snpres)
                    if (ns .eq. 0) then
                        call rctres(simec, snmec)
                    else
                        call rc32s0b(zr(jsnseis), simec, snmec)
                    end if
                    instp(1) = zr(jtranp+50*(inst1-1)+1)
                    instp(2) = -1.d0
!
                end if
            else
! ----- la première situation n'a pas de transitoire
                if (tranq) then
! ------- la deuxième situation a un transitoire
                    do l = 1, nbinstq
                        do i0 = 1, 2
                            do j = 1, 6
                                sc(j) = 0.d0
                                sb(j) = 0.d0-zr(jtranq+50*(l-1)+1+12+j)
                                sbet(j) = -( &
                                          zr(jtranq+50*(l-1)+1+12+j)-zr(jtranq+50*(l-1)+1+42+j))
                                do k = 1, 12
                                    if (nmecaq .eq. 3) mom2(k) = A1q(k)*zr(jtranq+50*(l-1)+1+48+1&
                                                                 &)+B1q(k)
                                    mij(k) = momcstp(k)-mom2(k)
                                    if (nmecap .eq. 3 .or. nmecaq .eq. 3) then
                                        sc(j) = sc(j)+mij(k)*zr(jsigu-1+78+6*(k-1)+j)
                                    end if
                                end do
                                st1(j) = sb(j)+sc(j)+e0(i0)*sa1(j)
                                stet1(j) = sbet(j)+sc(j)+e0(i0)*sa1(j)
                                if (unitaire) then
                                    st2(j) = sb(j)+sc(j)+e0(i0)*sa2(j)
                                    st3(j) = sb(j)+sc(j)+e0(i0)*sa3(j)
                                    st4(j) = sb(j)+sc(j)+e0(i0)*sa4(j)
                                    stet2(j) = sbet(j)+sc(j)+e0(i0)*sa2(j)
                                    stet3(j) = sbet(j)+sc(j)+e0(i0)*sa3(j)
                                    stet4(j) = sbet(j)+sc(j)+e0(i0)*sa4(j)
                                end if
                            end do
!
!
!------------------ Calcul du SN sans séisme
                            if (ns .eq. 0) then
                                call rctres(st1, tresca1)
                                if (unitaire) then
                                    call rctres(st2, tresca2)
                                    call rctres(st3, tresca3)
                                    call rctres(st4, tresca4)
                                end if
                                tresca = max(tresca1, tresca2, tresca3, tresca4)
!
                                if (tresca .gt. trescamax) then
                                    trescamax = tresca
                                    inst2 = l
                                end if
                            end if
!------------------ Calcul du SN avec séisme
                            if (ns .ne. 0) then
                                call rc32s0b(zr(jsnseis), st1, tresca1)
                                if (unitaire) then
                                    call rc32s0b(zr(jsnseis), st2, tresca2)
                                    call rc32s0b(zr(jsnseis), st3, tresca3)
                                    call rc32s0b(zr(jsnseis), st4, tresca4)
                                end if
                                tresca = max(tresca1, tresca2, tresca3, tresca4)
!
                                if (tresca .gt. trescamax) then
                                    trescamax = tresca
                                    inst2 = l
                                end if
                            end if
!------------------ Calcul du SN* sans séisme
                            if (ns .eq. 0) then
                                call rctres(stet1, tresca1)
                                if (unitaire) then
                                    call rctres(stet2, tresca2)
                                    call rctres(stet3, tresca3)
                                    call rctres(stet4, tresca4)
                                end if
                                tresca = max(tresca1, tresca2, tresca3, tresca4)
                                if (tresca .gt. trescaetmax) then
                                    trescaetmax = tresca
                                    instp(3) = -1.d0
                                    instp(4) = zr(jtranq+50*(l-1)+1)
                                end if
                            end if
!------------------ Calcul du SN* avec séisme
                            if (ns .ne. 0) then
                                call rc32s0b(zr(jsnseis), stet1, tresca1)
                                if (unitaire) then
                                    call rc32s0b(zr(jsnseis), stet2, tresca2)
                                    call rc32s0b(zr(jsnseis), stet3, tresca3)
                                    call rc32s0b(zr(jsnseis), stet4, tresca4)
                                end if
                                tresca = max(tresca1, tresca2, tresca3, tresca4)
                                if (tresca .gt. trescaetmax) then
                                    trescaetmax = tresca
                                    instp(3) = -1.d0
                                    instp(4) = zr(jtranq+50*(l-1)+1)
                                end if
                            end if
                            if (.not. unitaire) exit
!
                        end do
                    end do
!
!---------- Calcul de snther, snpres et snmec
                    do j = 1, 6
                        sith(j) = -zr(jtranq+50*(inst2-1)+1+24+j)
                        sipr(j) = -zr(jtranq+50*(inst2-1)+1+30+j)
                        simec(j) = -zr(jtranq+50*(inst2-1)+1+36+j)
                    end do
                    call rctres(sith, snther)
                    call rctres(sipr, snpres)
                    if (ns .eq. 0) then
                        call rctres(simec, snmec)
                    else
                        call rc32s0b(zr(jsnseis), simec, snmec)
                    end if
                    instp(1) = -1.d0
                    instp(2) = zr(jtranq+50*(inst2-1)+1)
!
                else
                    snther = 0.d0
                    snpres = 0.d0
                    snmec = 0.d0
                    instp(1) = -1.d0
                    instp(2) = -1.d0
                    instp(3) = -1.d0
                    instp(4) = -1.d0
! ------- la deuxième situation n'a pas non plus de transitoire
                    do i0 = 1, 2
                        do j = 1, 6
                            sc(j) = 0.d0
                            do k = 1, 12
                                mij(k) = momcstp(k)-momcstq(k)
                                if (nmecap .eq. 3 .or. nmecaq .eq. 3) then
                                    sc(j) = sc(j)+mij(k)*zr(jsigu-1+78+6*(k-1)+j)
                                end if
                            end do
                            st1(j) = sc(j)+e0(i0)*sa1(j)
                            st2(j) = sc(j)+e0(i0)*sa2(j)
                            st3(j) = sc(j)+e0(i0)*sa3(j)
                            st4(j) = sc(j)+e0(i0)*sa4(j)
                        end do
!
!
!------------------ Calcul du SN sans et avec séisme
                        if (ns .eq. 0) then
                            call rctres(st1, tresca1)
                            call rctres(st2, tresca2)
                            call rctres(st3, tresca3)
                            call rctres(st4, tresca4)
                            tresca = max(tresca1, tresca2, tresca3, tresca4)
                        else
                            call rc32s0b(zr(jsnseis), st1, tresca1)
                            call rc32s0b(zr(jsnseis), st2, tresca2)
                            call rc32s0b(zr(jsnseis), st3, tresca3)
                            call rc32s0b(zr(jsnseis), st4, tresca4)
                            tresca = max(tresca1, tresca2, tresca3, tresca4)
                        end if
!
                        if (tresca .gt. trescamax) then
                            trescamax = tresca
                            trescaetmax = tresca
                        end if
!
                    end do
                end if
            end if
!
            sn = trescamax
            snet = trescaetmax
!
        end if
!
        instsn(1) = instp(1)
        instsn(2) = instp(2)
        instsn(3) = instp(3)
        instsn(4) = instp(4)
    end if
!
!-----------------------------------------------------------
!   CALCUL DE LA PARTIE SP3 (QUI DEPEND DES INSTANTS DE SN)
!   POUR LE TERME DE RAPPEL AVEC LES INDICES DE CONTRAINTES
!-----------------------------------------------------------
!
    call jeveuo('&&RC3200.INDI', 'L', jvalin)
    k1 = zr(jvalin)
    c1 = zr(jvalin+1)
    k2 = zr(jvalin+2)
    c2 = zr(jvalin+3)
    k3 = zr(jvalin+4)
    c3 = zr(jvalin+5)
!
    sp3 = (k3*c3-1)*snther+(k1*c1-1)*snpres+(k2*c2-1)*snmec
    spmeca3 = (k1*c1-1)*snpres+(k2*c2-1)*snmec
!
!
    call jedema()
!
end subroutine
