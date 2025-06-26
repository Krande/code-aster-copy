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
subroutine rc32ac(lfat, lefat)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/r8vide.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rc32fact.h"
#include "asterfort/rc32pmb.h"
#include "asterfort/rc32s0.h"
#include "asterfort/rc32sa.h"
#include "asterfort/rc32sn.h"
#include "asterfort/rc32sp.h"
#include "asterfort/wkvect.h"
!
    aster_logical :: lfat, lefat
!
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE B3200 et ZE200
!                           CALCUL DES GRANDEURS
!
!     ------------------------------------------------------------------
    integer(kind=8) :: nb, ndim, jresu, iocc, im, jmax, ns, jresus, n1, i, jcombi
    integer(kind=8) :: iocc1, iocc2, jresucomb, jresucombs, jvalin, nbsscyc, nbid
    integer(kind=8) :: kk, jsnseis, jspseis, jinfo, numsit1, numsit2
    real(kind=8) :: pm, pb, pmpb, pmmax, pbmax, pmpbmax, pms, pbs, pmpbs
    real(kind=8) :: snmax, sn, sns, instsn(4), instsns(4), snet, snets
    real(kind=8) :: snetmax, sigmoypres, siprmoymax, snthmax, snther, sbid
    real(kind=8) :: spmax, sp(2), sps(2), instsp(4), instsps(4), samax
    real(kind=8) :: spmeca(2), spmecas(2), salt(2), salts(2), kemeca
    real(kind=8) :: kether, fu(2), fus(2), snp, snq, snetp, snetq, sigmoypresp
    real(kind=8) :: sigmoypresq, sntherp, sntherq, snthers, snps, snqs, snetps
    real(kind=8) :: snetqs, sntherps, sntherqs, spp, spq, futot, spthmax
    real(kind=8) :: fuseism, spmecap, spmecaq, futotenv, ktsn, ktsp, sp3
    real(kind=8) :: spmeca3, sp3s, spmeca3s, spss2(2), keequi, saltss(2)
    real(kind=8) :: fuss(2), spss(1000), spssbid(1000), fusstot
    character(len=4) :: lieu(2)
    character(len=16) :: option, chap
    character(len=8) :: typeke, sscyc
    aster_logical :: ze200
!
! DEB ------------------------------------------------------------------
!
    call jemarq()
!
    lfat = .false.
    lefat = .false.
    lieu(1) = 'ORIG'
    lieu(2) = 'EXTR'
!
!-- fait-on du calcul de pmpb, sn ou fatigue
    call getvtx(' ', 'OPTION', iocc=1, scal=option, nbret=n1)
!-- est on en ZE200
    ze200 = .false.
    call getvtx(' ', 'TYPE_RESU_MECA', iocc=1, scal=chap, nbret=n1)
    if (chap .eq. 'ZE200a' .or. chap .eq. 'ZE200b') ze200 = .true.
!-- est on en fatigue environnementale ?
    if (option .eq. 'EFAT') lefat = .true.
!-- combien de situations sont déclarées
    call getfac('SITUATION', nb)
!-- le séisme est-il pris en compte
    call getfac('SEISME', ns)
!
!-- repérer une situationpar son numéro dans le .mess
    call jeveuo('&&RC3200.SITU_INFOI', 'L', jinfo)
!
!-- les grandeurs sont calculées à l'origine et à l'extrémité du segment
    do im = 1, 2
!
! ---------------------------------------------------------------------------------
! - POUR TOUT LE CALCUL (SITUATION, COMBINAISON, SEISME) STOCKAGE DES 9 GRANDEURS
! - MAXIMALES A L'ORIGINE ET A L'EXTREMITE (&&RC3200.MAX_RESU)
! ---------------------------------------------------------------------------------
        pmmax = 0.d0
        pbmax = 0.d0
        pmpbmax = 0.d0
        snmax = 0.d0
        snetmax = 0.d0
        siprmoymax = 0.d0
        snthmax = 0.d0
        spmax = 0.d0
        samax = 0.d0
        spthmax = 0.d0
        fuseism = 0.d0
!
        call wkvect('&&RC3200.MAX_RESU.'//lieu(im), 'V V R', 13, jmax)
        do i = 1, 13
            zr(jmax+i-1) = r8vide()
        end do
! ---------------------------------------------------------------------------------
! - POUR CHAQUE SITUATION, STOCKAGE DE 123 GRANDEURS A L'ORIGINE ET A L'EXTREMITE
! - (&&RC3200.SITU_RESU SANS SEISME et &&RC3200.SITUS_RESU AVEC SEISME)
! - PM, PB, PMPB, SN, INST_SN_1, INST_SN_2, SN*, INST_SN*_1, INST_SN*_2
! - SIGMOYPRES, SN_THER, SP1, INST_SALT1_1, INST_SALT1_2, SALT, FUEL, SP_THER, SPMECA,
! - KE POUR EFAT, SPSS(100), NBSSCYC, FUSS, KE_MECA, KE_THER
! ---------------------------------------------------------------------------------
        ndim = nb*123
        call wkvect('&&RC3200.SITU_RESU.'//lieu(im), 'V V R', ndim, jresu)
        do i = 1, ndim
            zr(jresu+i-1) = r8vide()
        end do
        if (ns .ne. 0) then
            call wkvect('&&RC3200.SITUS_RESU.'//lieu(im), 'V V R', ndim, jresus)
            do i = 1, ndim
                zr(jresus+i-1) = r8vide()
            end do
        end if
!
!---- Calcul du séisme (contraintes linéarisées et totales) une seule fois
!
        if (ns .ne. 0 .and. .not. ze200) then
            call wkvect('&&RC3200.SNSEISME.'//lieu(im), 'V V R', 72, jsnseis)
            call rc32s0('SNSN', lieu(im), zr(jsnseis))
        end if
        if (option .eq. 'FATIGUE' .or. option .eq. 'EFAT') then
            if (ns .ne. 0 .and. .not. ze200) then
                call wkvect('&&RC3200.SPSEISME.'//lieu(im), 'V V R', 72, jspseis)
                call rc32s0('SPSP', lieu(im), zr(jspseis))
            end if
        end if
!
        do iocc = 1, nb, 1
!
            numsit1 = zi(jinfo+27*(iocc-1))
            if (im .eq. 1) write (*, 210) numsit1
            if (im .eq. 2) write (*, 211) numsit1
!
!------ Calcul du PM, PB et PMPB
            if (option .eq. 'PM_PB') then
                pm = 0.d0
                pb = 0.d0
                pmpb = 0.d0
                pms = 0.d0
                pbs = 0.d0
                pmpbs = 0.d0
!
                call rc32pmb(lieu(im), iocc, 0, pm, pb, &
                             pmpb)
!
                zr(jresu+123*(iocc-1)) = pm
                zr(jresu+123*(iocc-1)+1) = pb
                zr(jresu+123*(iocc-1)+2) = pmpb
                pmmax = max(pm, pmmax)
                pbmax = max(pb, pbmax)
                pmpbmax = max(pmpb, pmpbmax)
!
                if (ns .ne. 0) then
                    call rc32pmb(lieu(im), iocc, ns, pms, pbs, &
                                 pmpbs)
                    zr(jresus+123*(iocc-1)) = pms
                    zr(jresus+123*(iocc-1)+1) = pbs
                    zr(jresus+123*(iocc-1)+2) = pmpbs
                    pmmax = max(pms, pmmax)
                    pbmax = max(pbs, pbmax)
                    pmpbmax = max(pmpbs, pmpbmax)
                end if
            end if
!
!------ Calcul du SN
            if (option .eq. 'SN' .or. option .eq. 'FATIGUE' .or. option .eq. 'EFAT') then
                call jeveuo('&&RC3200.INDI', 'L', jvalin)
                ktsn = zr(jvalin+9)
                ktsp = zr(jvalin+10)
                sn = 0.d0
                sns = 0.d0
                snet = 0.d0
                snets = 0.d0
                snther = 0.d0
                sp3 = 0.d0
                spmeca3 = 0.d0
                sp3s = 0.d0
                spmeca3s = 0.d0
!
                call rc32sn(ze200, lieu(im), iocc, 0, 0, &
                            sn, instsn, snet, sigmoypres, snther, &
                            sp3, spmeca3)
                sn = sn*ktsn
                write (*, *) 'SN=', sn
                snet = snet*ktsn
                zr(jresu+123*(iocc-1)+3) = sn
                zr(jresu+123*(iocc-1)+4) = instsn(1)
                zr(jresu+123*(iocc-1)+5) = instsn(2)
                zr(jresu+123*(iocc-1)+6) = snet
                zr(jresu+123*(iocc-1)+7) = instsn(3)
                zr(jresu+123*(iocc-1)+8) = instsn(4)
                zr(jresu+123*(iocc-1)+9) = sigmoypres
                zr(jresu+123*(iocc-1)+10) = snther
                snmax = max(sn, snmax)
                snetmax = max(snet, snetmax)
                snthmax = max(snther, snthmax)
                if (sigmoypres .ne. r8vide()) siprmoymax = max(sigmoypres, siprmoymax)
!
                if (ns .ne. 0) then
                    call rc32sn(ze200, lieu(im), iocc, 0, ns, &
                                sns, instsns, snets, sbid, sbid, &
                                sp3s, spmeca3s)
                    sns = ktsn*sns
                    write (*, *) 'SN=', sns
                    snets = ktsn*snets
                    zr(jresus+123*(iocc-1)+3) = sns
                    zr(jresus+123*(iocc-1)+4) = instsns(1)
                    zr(jresus+123*(iocc-1)+5) = instsns(2)
                    zr(jresus+123*(iocc-1)+6) = snets
                    zr(jresus+123*(iocc-1)+7) = instsns(3)
                    zr(jresus+123*(iocc-1)+8) = instsns(4)
                    zr(jresus+123*(iocc-1)+9) = sigmoypres
                    zr(jresus+123*(iocc-1)+10) = snther
                    snmax = max(sns, snmax)
                    snetmax = max(snets, snetmax)
                end if
            end if
!
!------ Calcul du SP, du SALT et du FU de chaque situation seule
            if (option .eq. 'FATIGUE' .or. option .eq. 'EFAT') then
                lfat = .true.
                sp(1) = 0.d0
                salt(1) = 0.d0
                spmeca(1) = 0.d0
                sps(1) = 0.d0
                salts(1) = 0.d0
                spmecas(1) = 0.d0
                instsp(1) = -1.d0
                instsp(2) = -1.d0
                instsps(1) = -1.d0
                instsps(2) = -1.d0
!
                call rc32sp(ze200, lieu(im), iocc, 0, 0, &
                            sp, spmeca, instsp, nbsscyc, spss)
                sp(1) = ktsp*(sp(1)+sp3)
                write (*, *) 'SP=', sp(1)
                spmeca(1) = ktsp*(spmeca(1)+spmeca3)
                zr(jresu+123*(iocc-1)+11) = sp(1)
                zr(jresu+123*(iocc-1)+12) = instsp(1)
                zr(jresu+123*(iocc-1)+13) = instsp(2)
                spmax = max(sp(1), spmax)
                call rc32sa('SITU', sn, sp, spmeca, kemeca, &
                            kether, salt, fu)
                zr(jresu+123*(iocc-1)+14) = salt(1)
                zr(jresu+123*(iocc-1)+15) = fu(1)
                zr(jresu+123*(iocc-1)+16) = max(0.d0, sp(1)-spmeca(1))
                zr(jresu+123*(iocc-1)+17) = spmeca(1)
!
                call getvtx(' ', 'TYPE_KE', scal=typeke, nbret=n1)
                if (typeke .eq. 'KE_MECA') then
                    keequi = kemeca
                else
                    keequi = (kemeca*spmeca(1)+kether*(sp(1)-spmeca(1)))/sp(1)
                end if
                zr(jresu+123*(iocc-1)+18) = keequi
                zr(jresu+123*(iocc-1)+121) = kemeca
                zr(jresu+123*(iocc-1)+122) = kether
!
!---------- Calcul du FU dû aux sous cycles
                call getvtx(' ', 'SOUS_CYCL', scal=sscyc, nbret=n1)
                if (sscyc .eq. 'OUI') then
                    fusstot = 0.d0
                    do kk = 1, nbsscyc
                        zr(jresu+123*(iocc-1)+19-1+kk) = spss(kk)
                        spss2(1) = keequi*spss(kk)
                        spss2(2) = 0.d0
                        call rc32sa('SITU', 1.d-12, spss2, spss2, kemeca, &
                                    kether, saltss, fuss)
                        fusstot = fusstot+fuss(1)
                    end do
                    zr(jresu+123*(iocc-1)+119) = nbsscyc
                    zr(jresu+123*(iocc-1)+120) = fusstot
                    zr(jresu+123*(iocc-1)+15) = zr(jresu+123*(iocc-1)+15)+fusstot
                else
                    zr(jresu+123*(iocc-1)+120) = 0.d0
                end if
!
                samax = max(salt(1), samax)
                spthmax = max(spthmax, zr(jresu+123*(iocc-1)+16))
!
!------ Calcul du SP, du SALT et du FU avec SEISME
                if (ns .ne. 0) then
                    call rc32sp(ze200, lieu(im), iocc, 0, ns, &
                                sps, spmecas, instsps, nbid, spssbid)
                    sps(1) = ktsp*(sps(1)+sp3s)
                    write (*, *) 'SPS=', sps(1)
                    spmecas(1) = ktsp*(spmecas(1)+spmeca3s)
                    zr(jresus+123*(iocc-1)+11) = sps(1)
                    zr(jresus+123*(iocc-1)+12) = instsps(1)
                    zr(jresus+123*(iocc-1)+13) = instsps(2)
                    spmax = max(spmax, sps(1))
                    call rc32sa('SITU', sns, sps, spmecas, kemeca, &
                                kether, salts, fus)
                    zr(jresus+123*(iocc-1)+14) = salts(1)
                    zr(jresus+123*(iocc-1)+15) = fus(1)
                    zr(jresus+123*(iocc-1)+16) = max(0.d0, sps(1)-spmecas(1))
                    zr(jresus+123*(iocc-1)+17) = spmecas(1)
!
                    if (typeke .eq. 'KE_MECA') then
                        keequi = kemeca
                    else
                        keequi = (kemeca*spmecas(1)+kether*(sps(1)-spmecas(1)))/sps(1)
                    end if
                    zr(jresus+123*(iocc-1)+18) = keequi
                    zr(jresus+123*(iocc-1)+121) = kemeca
                    zr(jresus+123*(iocc-1)+122) = kether
!
!
!---------- Calcul du FU dû aux sous cycles avec SEISME
!
                    if (sscyc .eq. 'OUI') then
                        fusstot = 0.d0
                        do kk = 1, nbsscyc
                            zr(jresus+123*(iocc-1)+19-1+kk) = spss(kk)
                            spss2(1) = keequi*spss(kk)
                            spss2(2) = 0.d0
                            call rc32sa('SITU', 1.d-12, spss2, spss2, kemeca, &
                                        kether, saltss, fuss)
                            fusstot = fusstot+fuss(1)
                        end do
                        zr(jresus+123*(iocc-1)+119) = nbsscyc
                        zr(jresus+123*(iocc-1)+120) = fusstot
                        zr(jresus+123*(iocc-1)+15) = zr(jresu+123*(iocc-1)+15)+fusstot
                    else
                        zr(jresus+123*(iocc-1)+120) = 0.d0
                    end if
!
                    samax = max(samax, salts(1))
                    spthmax = max(spthmax, zr(jresus+123*(iocc-1)+16))
!
                end if
            end if
!
        end do
!
! ---------------------------------------------------------------------------------
! - POUR CHAQUE COMBINAISON DE SITUATIONS POSSIBLE (MEME GROUPE ou SITUATION DE PASSAGE)
! - STOCKAGE DE 25 GRANDEURS A L'ORIGINE ET A L'EXTREMITE
! - (&&RC3200.SITU_RESU SANS SEISME et &&RC3200.SITUS_RESU AVEC SEISME)
! - SN, INST_SN_1, INST_SN_2, SN*, INST_SN*_1, INST_SN*_2
! - SIGMOYPRES, SN_THER, SP1,INST_SALT1_1, INST_SALT1_2, SALT1, SP2,INST_SALT2_1,
! - INST_SALT2_2, SALT2, FUEL, SP_THER1(18), SP_THER2(19), KE POUR EFAT(20), FUSS(21)
! - , KE_MECA(22), KE_THER(23), FUEL1(24),FUEL2(25)
! ---------------------------------------------------------------------------------
!
        ndim = nb*nb*25
        call wkvect('&&RC3200.COMB_RESU.'//lieu(im), 'V V R', ndim, jresucomb)
        do i = 1, ndim
            zr(jresucomb+i-1) = r8vide()
        end do
        if (ns .ne. 0) then
            call wkvect('&&RC3200.COMBS_RESU.'//lieu(im), 'V V R', ndim, jresucombs)
            do i = 1, ndim
                zr(jresucombs+i-1) = r8vide()
            end do
        end if
!
!-- on consulte le tableau qui indique si la combinaison des deux situations P et Q
!-- est possible (elles sont dans le même groupe ou liées par une situation de passage)
        if (option .eq. 'FATIGUE' .or. option .eq. 'EFAT') then
            call jeveuo('&&RC3200.COMBI', 'L', jcombi)
            do iocc1 = 1, nb
                do iocc2 = 1, nb
                    sn = 0.d0
                    sns = 0.d0
                    snet = 0.d0
                    snets = 0.d0
                    snther = 0.d0
                    if (zi(jcombi+nb*(iocc1-1)+iocc2-1) .ne. 0) then
!
!---------------- Calcul de SN(P,Q), SN*(P,Q) et leurs instants sans séisme
                        numsit1 = zi(jinfo+27*(iocc1-1))
                        numsit2 = zi(jinfo+27*(iocc2-1))
                        if (im .eq. 1) write (*, 110) numsit1, numsit2
                        if (im .eq. 2) write (*, 111) numsit1, numsit2
                        call rc32sn(ze200, lieu(im), iocc1, iocc2, 0, &
                                    sn, instsn, snet, sbid, snther, &
                                    sp3, spmeca3)
                        sn = ktsn*sn
                        write (*, *) 'SN=', sn
                        snet = ktsn*snet
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+1) = sn
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+2) = instsn(1)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+3) = instsn(2)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+4) = snet
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+5) = instsn(3)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+6) = instsn(4)
                        sigmoypresp = zr(jresu+123*(iocc1-1)+9)
                        sigmoypresq = zr(jresu+123*(iocc2-1)+9)
                      zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+7) = max(sigmoypresp, sigmoypresq)
                        if (sigmoypresp .ne. r8vide()) siprmoymax = max(sigmoypresp, siprmoymax)
                        if (sigmoypresq .ne. r8vide()) siprmoymax = max(sigmoypresq, siprmoymax)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+8) = snther
!---------------- On prend le max de SN(P,Q), SN(P,P) et SN(Q,Q)
                        snp = zr(jresu+123*(iocc1-1)+3)
                        snq = zr(jresu+123*(iocc2-1)+3)
                        sntherp = zr(jresu+123*(iocc1-1)+10)
                        sntherq = zr(jresu+123*(iocc2-1)+10)
                        if (snp .ge. sn) then
                            sn = snp
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+1) = sn
                          zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+2) = zr(jresu+123*(iocc1-1)+4)
                          zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+3) = zr(jresu+123*(iocc1-1)+5)
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+8) = sntherp
                        end if
                        if (snq .ge. sn) then
                            sn = snq
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+1) = sn
                          zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+2) = zr(jresu+123*(iocc2-1)+4)
                          zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+3) = zr(jresu+123*(iocc2-1)+5)
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+8) = sntherq
                        end if
                        snmax = max(snmax, sn)
                        snthmax = max(zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+8), snthmax)
!
                        snetp = zr(jresu+123*(iocc1-1)+6)
                        snetq = zr(jresu+123*(iocc2-1)+6)
                        if (snetp .ge. snet) then
                            snet = snetp
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+4) = snet
                          zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+5) = zr(jresu+123*(iocc1-1)+7)
                          zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+6) = zr(jresu+123*(iocc1-1)+8)
                        end if
                        if (snetq .ge. snet) then
                            snet = snetq
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+4) = snet
                          zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+5) = zr(jresu+123*(iocc2-1)+7)
                          zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+6) = zr(jresu+123*(iocc2-1)+8)
                        end if
                        snetmax = max(snet, snetmax)
!
!---------------- Calcul de SN, SN* et leurs instants avec séisme
                        if (ns .ne. 0) then
                            call rc32sn(ze200, lieu(im), iocc1, iocc2, ns, &
                                        sns, instsns, snets, sbid, snthers, &
                                        sp3s, spmeca3s)
                            sns = ktsn*sns
                            write (*, *) 'SNS=', sns
                            snets = ktsn*snets
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+1) = sns
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+2) = instsns(1)
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+3) = instsns(2)
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+4) = snets
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+5) = instsns(3)
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+6) = instsns(4)
                     zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+7) = max(sigmoypresp, sigmoypresq)
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+8) = snthers
!
                            snps = zr(jresus+123*(iocc1-1)+3)
                            snqs = zr(jresus+123*(iocc2-1)+3)
                            sntherps = zr(jresus+123*(iocc1-1)+10)
                            sntherqs = zr(jresus+123*(iocc2-1)+10)
                            if (snps .ge. sns) then
                                sns = snps
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+1) = sns
                        zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+2) = zr(jresus+123*(iocc1-1)+4)
                        zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+3) = zr(jresus+123*(iocc1-1)+5)
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+8) = sntherps
                            end if
                            if (snqs .ge. sns) then
                                sns = snqs
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+1) = sns
                        zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+2) = zr(jresus+123*(iocc2-1)+4)
                        zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+3) = zr(jresus+123*(iocc2-1)+5)
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+8) = sntherqs
                            end if
                            snmax = max(snmax, sns)
                            snthmax = max(zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+8), snthmax)
!
                            snetps = zr(jresus+123*(iocc1-1)+6)
                            snetqs = zr(jresus+123*(iocc2-1)+6)
                            if (snetps .ge. snets) then
                                snets = snetps
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+4) = snets
                        zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+5) = zr(jresus+123*(iocc1-1)+7)
                        zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+6) = zr(jresus+123*(iocc1-1)+8)
                            end if
                            if (snetqs .ge. snets) then
                                snets = snetqs
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+4) = snets
                        zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+5) = zr(jresus+123*(iocc2-1)+7)
                        zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+6) = zr(jresus+123*(iocc2-1)+8)
                            end if
                            snetmax = max(snets, snetmax)
                        end if
!
!---------------- Calcul de SP1(P,Q), SP2(P,Q) et leurs instants sans séisme
                        call rc32sp(ze200, lieu(im), iocc1, iocc2, 0, &
                                    sp, spmeca, instsp, nbid, spssbid)
                        sp(1) = ktsp*(sp(1)+sp3)
                        write (*, *) 'SP=', sp(1)
                        spmeca(1) = ktsp*(spmeca(1)+spmeca3)
                        sp(2) = ktsp*(sp(2)+sp3)
                        spmeca(2) = ktsp*(spmeca(2)+spmeca3)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+9) = sp(1)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+10) = instsp(1)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+11) = instsp(2)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+13) = sp(2)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+14) = instsp(3)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+15) = instsp(4)
                        call rc32sa('COMB', sn, sp, spmeca, kemeca, &
                                    kether, salt, fu)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+12) = salt(1)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+16) = salt(2)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+17) = fu(1)+fu(2)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+18) = max(0.d0, sp(1)-spmeca(1))
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+19) = max(0.d0, sp(2)-spmeca(2))
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+24) = fu(1)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+25) = fu(2)
!
!---------------- On prend le max de SP(P,Q), SP(P,P) et SP(Q,Q)
                        spp = zr(jresu+123*(iocc1-1)+11)
                        spq = zr(jresu+123*(iocc2-1)+11)
                        spmecap = zr(jresu+123*(iocc1-1)+17)
                        spmecaq = zr(jresu+123*(iocc2-1)+17)
                        if (spp .ge. sp(1)) then
                            sp(1) = spp
                            sp(2) = spq
                            spmeca(1) = spmecap
                            spmeca(2) = spmecaq
                            call rc32sa('COMB', sn, sp, spmeca, kemeca, &
                                        kether, salt, fu)
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+9) = spp
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+10) = zr(jresu+123*(iocc1-1)+12)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+11) = zr(jresu+123*(iocc1-1)+13)
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+12) = salt(1)
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+13) = spq
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+14) = zr(jresu+123*(iocc2-1)+12)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+15) = zr(jresu+123*(iocc2-1)+13)
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+16) = salt(2)
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+17) = fu(1)+fu(2)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+18) = max(0.d0, sp(1)-spmeca(1))
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+19) = max(0.d0, sp(2)-spmeca(2))
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+24) = fu(1)
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+25) = fu(2)
                        end if
                        if (spq .ge. sp(1)) then
                            sp(1) = spq
                            sp(2) = spp
                            spmeca(1) = spmecaq
                            spmeca(2) = spmecap
                            call rc32sa('COMB', sn, sp, spmeca, kemeca, &
                                        kether, salt, fu)
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+9) = spq
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+10) = zr(jresu+123*(iocc2-1)+12)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+11) = zr(jresu+123*(iocc2-1)+13)
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+12) = salt(1)
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+13) = spp
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+14) = zr(jresu+123*(iocc1-1)+12)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+15) = zr(jresu+123*(iocc1-1)+13)
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+16) = salt(2)
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+17) = fu(1)+fu(2)
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+18) = max(0.d0, sp(1)-spmeca(1))
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+19) = max(0.d0, sp(2)-spmeca(2))
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+24) = fu(1)
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+25) = fu(2)
                        end if
!
                        if (typeke .eq. 'KE_MECA') then
                            keequi = kemeca
                        else
                            keequi = (kemeca*spmeca(1)+kether*(sp(1)-spmeca(1)))/sp(1)
                        end if
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+20) = keequi
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+22) = kemeca
                        zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+23) = kether
!
!
!---------- Calcul du FU dû aux sous cycles
!
                        if (sscyc .eq. 'OUI') then
                            fusstot = 0.d0
                            nbsscyc = int(zr(jresu+123*(iocc1-1)+119))
                            do kk = 1, nbsscyc
                                spss2(1) = keequi*zr(jresu+123*(iocc1-1)+19-1+kk)
                                spss2(2) = 0.d0
                                call rc32sa('SITU', 1.d-12, spss2, spss2, kemeca, &
                                            kether, saltss, fuss)
                                fusstot = fusstot+fuss(1)
                            end do
                            nbsscyc = int(zr(jresu+123*(iocc2-1)+119))
                            do kk = 1, nbsscyc
                                spss2(1) = keequi*zr(jresu+123*(iocc2-1)+19-1+kk)
                                spss2(2) = 0.d0
                                call rc32sa('SITU', 1.d-12, spss2, spss2, kemeca, &
                                            kether, saltss, fuss)
                                fusstot = fusstot+fuss(1)
                            end do
!
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+21) = fusstot
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+17) = &
                                zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+17)+fusstot
                        else
                            zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+21) = 0.d0
                        end if
!
                        spmax = max(spmax, sp(1))
                        samax = max(samax, salt(1))
                        spthmax = max(spthmax, zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+18))
!
!---------------- Calcul de SP1(P,Q), SP2(P,Q) et leurs instants avec séisme
                        if (ns .ne. 0) then
                            call rc32sp(ze200, lieu(im), iocc1, iocc2, ns, &
                                        sps, spmecas, instsps, nbid, spssbid)
                            sps(1) = ktsp*(sps(1)+sp3s)
                            write (*, *) 'SPS=', sps(1)
                            spmecas(1) = ktsp*(spmecas(1)+spmeca3s)
                            sps(2) = ktsp*(sps(2)+sp3s)
                            spmecas(2) = ktsp*(spmecas(2)+spmeca3s)
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+9) = sps(1)
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+10) = instsps(1)
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+11) = instsps(2)
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+13) = sps(2)
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+14) = instsps(3)
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+15) = instsps(4)
                            call rc32sa('COMB', sns, sps, spmecas, kemeca, &
                                        kether, salts, fus)
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+12) = salts(1)
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+16) = salts(2)
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+17) = fus(1)+fus(2)
                     zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+18) = max(0.d0, sps(1)-spmecas(1))
                     zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+19) = max(0.d0, sps(2)-spmecas(2))
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+24) = fus(1)
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+25) = fus(2)
!---------------- On prend le max de SP(P,Q), SP(P,P) et SP(Q,Q)
                            spp = zr(jresus+123*(iocc1-1)+11)
                            spq = zr(jresus+123*(iocc2-1)+11)
                            spmecap = zr(jresus+123*(iocc1-1)+17)
                            spmecaq = zr(jresus+123*(iocc2-1)+17)
                            if (spp .ge. sps(1)) then
                                sps(1) = spp
                                sps(2) = spq
                                spmecas(1) = spmecap
                                spmecas(2) = spmecaq
                                call rc32sa('COMB', sns, sps, spmecas, kemeca, &
                                            kether, salts, fus)
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+9) = spp
                      zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+10) = zr(jresus+123*(iocc1-1)+12)
                      zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+11) = zr(jresus+123*(iocc1-1)+13)
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+12) = salts(1)
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+13) = spq
                      zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+14) = zr(jresus+123*(iocc2-1)+12)
                      zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+15) = zr(jresus+123*(iocc2-1)+13)
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+16) = salts(2)
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+17) = fus(1)+fus(2)
                     zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+18) = max(0.d0, sps(1)-spmecas(1))
                     zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+19) = max(0.d0, sps(2)-spmecas(2))
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+24) = fus(1)
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+25) = fus(2)
                            end if
                            if (spq .ge. sps(1)) then
                                sps(1) = spq
                                sps(2) = spp
                                spmecas(1) = spmecaq
                                spmecas(2) = spmecap
                                call rc32sa('COMB', sns, sps, spmecas, kemeca, &
                                            kether, salts, fus)
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+9) = spq
                      zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+10) = zr(jresus+123*(iocc2-1)+12)
                      zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+11) = zr(jresus+123*(iocc2-1)+13)
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+12) = salts(1)
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+13) = spp
                      zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+14) = zr(jresus+123*(iocc1-1)+12)
                      zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+15) = zr(jresus+123*(iocc1-1)+13)
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+16) = salts(2)
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+17) = fus(1)+fus(2)
                     zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+18) = max(0.d0, sps(1)-spmecas(1))
                     zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+19) = max(0.d0, sps(2)-spmecas(2))
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+24) = fus(1)
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+25) = fus(2)
                            end if
!
                            if (typeke .eq. 'KE_MECA') then
                                keequi = kemeca
                            else
                                keequi = (kemeca*spmecas(1)+kether*(sps(1)-spmecas(1)))/sps(1)
                            end if
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+20) = keequi
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+22) = kemeca
                            zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+23) = kether
!
!
!---------- Calcul du FU dû aux sous cycles
!
                            if (sscyc .eq. 'OUI') then
                                fusstot = 0.d0
                                nbsscyc = int(zr(jresu+123*(iocc1-1)+119))
                                do kk = 1, nbsscyc
                                    spss2(1) = keequi*zr(jresu+123*(iocc1-1)+19-1+kk)
                                    spss2(2) = 0.d0
                                    call rc32sa('SITU', 1.d-12, spss2, spss2, kemeca, &
                                                kether, saltss, fuss)
                                    fusstot = fusstot+fuss(1)
                                end do
                                nbsscyc = int(zr(jresu+123*(iocc2-1)+119))
                                do kk = 1, nbsscyc
                                    spss2(1) = keequi*zr(jresu+123*(iocc2-1)+19-1+kk)
                                    spss2(2) = 0.d0
                                    call rc32sa('SITU', 1.d-12, spss2, spss2, kemeca, &
                                                kether, saltss, fuss)
                                    fusstot = fusstot+fuss(1)
                                end do
!
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+21) = fusstot
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+17) = &
                                    zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+17)+fusstot
                            else
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+21) = 0.d0
                            end if
!
                            spmax = max(spmax, sps(1))
                            samax = max(samax, salts(1))
                            spthmax = max( &
                                      spthmax, zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+18))
                        end if
!
                    end if
                end do
            end do
!
! ---------------------------------------------------------------------------------
!                          CALCUL DU FACTEUR D'USAGE TOTAL
! ---------------------------------------------------------------------------------
            call rc32fact(ze200, nb, lieu(im), ns, fuseism, &
                          futot, lefat, futotenv)
        end if
! ---------------------------------------------------------------------------------
! --- on remplit le vecteur des quantités max (&&RC3200.MAX_RESU)
! ---------------------------------------------------------------------------------
        if (option .eq. 'PM_PB') then
            zr(jmax) = pmmax
            zr(jmax+1) = pbmax
            zr(jmax+2) = pmpbmax
        end if
        if (option .eq. 'SN' .or. option .eq. 'FATIGUE' .or. option .eq. 'EFAT') then
            zr(jmax+3) = snmax
            zr(jmax+4) = snetmax
            zr(jmax+5) = siprmoymax
            zr(jmax+6) = snthmax
        end if
        if (option .eq. 'FATIGUE' .or. option .eq. 'EFAT') then
            zr(jmax+7) = spthmax
            zr(jmax+8) = spmax
            zr(jmax+9) = samax
            zr(jmax+10) = futot
            zr(jmax+11) = fuseism
            zr(jmax+12) = futotenv
        end if
!
    end do
!
110 format(1p, ' ORIGINE, COMBINAISON DES SITUATIONS', i4, 3x, i4)
111 format(1p, ' EXTREMITE, COMBINAISON DES SITUATIONS', i4, 3x, i4)
210 format(1p, ' ORIGINE, SITUATION NUMERO', i4)
211 format(1p, ' EXTREMITE, SITUATION NUMERO', i4)
!
    call jedema()
!
end subroutine
