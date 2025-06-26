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
subroutine rc32rs(lfat, lefat)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterc/r8vide.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rcmcrt.h"
#include "asterfort/rcvale.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
!
    aster_logical :: lfat, lefat
!
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE B3200 et ZE200
!     AFFICHAGE DES RESULTATS DANS LA TABLE DE SORTIE
!
!     ------------------------------------------------------------------
    character(len=4) :: lieu(2)
    character(len=8) :: nomres, k8b, mater
    character(len=16) :: concep, nomcmd, typtab, valek(5)
    integer(kind=8) :: ibid, n1, npar1, npar0, npar2, im, valei(7), nb, i, jval
    integer(kind=8) :: ns, jnom, jresu, jmax, jresus, k, icodre(1), iocc1, iocc2
    integer(kind=8) :: jcombi, jresucomb, jresucombs, npar3, jfact, npar4
    integer(kind=8) :: num1, num2, noccpris, ind1, ind2, n5, jfactenv
    complex(kind=8) :: c16b
    parameter(npar0=55, npar1=18, npar2=29, npar3=29, npar4=32)
    character(len=16) :: nopar1(npar1), nopar0(npar0), nopar2(npar2)
    character(len=8) :: sscyc, typar1(npar1), typar0(npar0), typar4(npar4)
    character(len=16) :: nopar3(npar3), nopar4(npar4)
    real(kind=8) :: valer(23), sigmoypres, symax, rbid, valres(1)
    real(kind=8) :: fuseism, fenint, fenglobal
!
!  ------------------------------------------------------------------
    data lieu/'ORIG', 'EXTR'/
!
    data nopar0/'TYPE', 'SEISME', 'LIEU', 'NOM_SIT1',&
     &            'NUM_SIT1', 'NOCC_SIT1', 'GROUP_SIT1', 'PM', 'PB', 'PMPB',&
     &            'NOM_SIT2', 'NUM_SIT2', 'NOCC_SIT2', 'GROUP_SIT2', 'SN', 'INST_SN_1',&
     &            'INST_SN_2', 'SN*', 'INST_SN*_1', 'INST_SN*_2', &
     &            'SIG_PRES_MOY', 'SN_THER', 'CRIT_LINE', 'CRIT_PARAB',&
     &            'KE_MECA', 'KE_THER', 'SP(MECA)', 'SP(THER)', 'SP', &
     &            'INST_SALT_1', 'INST_SALT_2', 'SALT', 'FU_UNIT_SS_CYCL', &
     &            'FU_UNIT', 'NOCC_PRIS', 'FU_PARTIEL', 'FEN', 'FEN_ELAS',&
     &            'FUEN_PARTIEL', 'PM_MAX', 'PB_MAX', 'PMPB_MAX', 'SN_MAX',&
     &            'SN*_MAX', 'SIGM_M_PRES', 'SN_THER_MAX', 'CRIT_LINE_MAX',&
     &            'SP_THER_MAX', 'CRIT_PARA_MAX', 'SP_MAX', 'SALT_MAX', 'FU_TOTAL',&
     &            'FEN_TOTAL', 'FEN_INTEGRE', 'FUEN_TOTAL'/
    data typar0/'K8', 'K8', 'K8', 'K16',&
     &             'I', 'I', 'I', 'R', 'R', 'R',&
     &             'K16', 'I', 'I', 'I', 'R', 'R',&
     &             'R', 'R', 'R', 'R',&
     &             'R', 'R', 'R', 'R', &
     &             'R', 'R', 'R', 'R', 'R',&
     &             'R', 'R', 'R', 'R', &
     &             'R', 'I', 'R', 'R', 'R',&
     &             'R', 'R', 'R', 'R', 'R',&
     &             'R', 'R', 'R', 'R', &
     &             'R', 'R', 'R', 'R', 'R',&
     &             'R', 'R', 'R'/
!
    data nopar1/'TYPE', 'LIEU', 'PM_MAX', 'PB_MAX', 'PMPB_MAX', 'SN_MAX',&
     &            'SN*_MAX', 'SIGM_M_PRES', 'SN_THER_MAX', 'CRIT_LINE_MAX',&
     &            'SP_THER_MAX', 'CRIT_PARA_MAX', 'SP_MAX', 'SALT_MAX', 'FU_TOTAL',&
     &            'FEN_TOTAL', 'FEN_INTEGRE', 'FUEN_TOTAL'/
    data typar1/'K8', 'K8', 'R', 'R', 'R', 'R',&
     &            'R', 'R', 'R', 'R',&
     &            'R', 'R', 'R', 'R', 'R',&
     &            'R', 'R', 'R'/
!
    data nopar2/'TYPE', 'SEISME', 'LIEU',&
     &            'NOM_SIT1', 'NUM_SIT1', 'GROUP_SIT1', 'PM', 'PB', 'PMPB',&
     &            'SN', 'INST_SN_1', 'INST_SN_2', 'SN*', 'INST_SN*_1',&
     &            'INST_SN*_2', 'SIG_PRES_MOY', 'SN_THER', 'CRIT_LINE',&
     &            'CRIT_PARAB', 'KE_MECA', 'KE_THER', 'SP(MECA)', 'SP(THER)', 'SP',&
     &            'INST_SALT_1', 'INST_SALT_2', 'SALT', 'FU_UNIT_SS_CYCL', 'FU_UNIT'/
!
    data nopar3/'TYPE', 'SEISME', 'LIEU',&
     &            'NOM_SIT1', 'NUM_SIT1', 'GROUP_SIT1',&
     &            'NOM_SIT2', 'NUM_SIT2', 'GROUP_SIT2', 'SN',&
     &            'INST_SN_1', 'INST_SN_2', 'SN*', 'INST_SN*_1',&
     &            'INST_SN*_2', 'SIG_PRES_MOY', 'SN_THER', 'CRIT_LINE',&
     &            'CRIT_PARAB', 'KE_MECA', 'KE_THER', 'SP(MECA)', 'SP(THER)', &
     &            'SP', 'INST_SALT_1', 'INST_SALT_2', 'SALT', 'FU_UNIT_SS_CYCL', 'FU_UNIT'/
!
    data nopar4/'TYPE', 'SEISME', 'LIEU',&
     &            'NOM_SIT1', 'NUM_SIT1', 'NOCC_SIT1', 'GROUP_SIT1',&
     &            'NOM_SIT2', 'NUM_SIT2', 'NOCC_SIT2', 'GROUP_SIT2', 'SN',&
     &            'INST_SN_1', 'INST_SN_2', 'SN*', 'INST_SN*_1',&
     &            'INST_SN*_2', 'KE_MECA', 'KE_THER', 'SP(MECA)', 'SP(THER)', 'SP',&
     &            'INST_SALT_1', 'INST_SALT_2', 'SALT', 'FU_UNIT_SS_CYCL', 'FU_UNIT', &
                   'NOCC_PRIS', 'FU_PARTIEL', 'FEN', 'FEN_ELAS', 'FUEN_PARTIEL'/
    data typar4/'K8', 'K8', 'K8',&
     &            'K16', 'I', 'I', 'I',&
     &            'K16', 'I', 'I', 'I', 'R',&
     &            'R', 'R', 'R', 'R',&
     &            'R', 'R', 'R', 'R', 'R', 'R', &
     &            'R', 'R', 'R', 'R', 'R', 'I', &
     &            'R', 'R', 'R', 'R'/
!
!
! DEB ------------------------------------------------------------------
!
    call jemarq()
!
    call getres(nomres, concep, nomcmd)
!
    call tbcrsd(nomres, 'G')
!
    ibid = 0
    rbid = 0.d0
    valer = 0.d0
    c16b = (0.d0, 0.d0)
    call getvtx(' ', 'TYPE_RESU', scal=typtab, nbret=n1)
    call getvtx(' ', 'SOUS_CYCL', scal=sscyc, nbret=n1)
    call getvid(' ', 'MATER', scal=mater, nbret=n1)
!
!
    if (typtab .eq. 'VALE_MAX') then
        call tbajpa(nomres, npar1, nopar1, typar1)
    else if (typtab .eq. 'SYSTUS') then
        call tbajpa(nomres, npar4, nopar4, typar4)
    else
        call tbajpa(nomres, npar0, nopar0, typar0)
    end if
!     ----------------------------------------------------------------
! --- AFFICHAGE DES MAXIMA DANS LA TABLE
!     ----------------------------------------------------------------
    valek(1) = 'MAXI'
    do im = 1, 2
!
        valek(2) = lieu(im)
        call jeveuo('&&RC3200.MAX_RESU.'//lieu(im), 'L', jmax)
!
        do k = 1, 7
            valer(k) = zr(jmax-1+k)
        end do
!
        sigmoypres = zr(jmax+5)
        if (sigmoypres .eq. r8vide() .or. abs(sigmoypres) .lt. 1e-10) then
            valer(7) = r8vide()
            valer(8) = r8vide()
            valer(10) = r8vide()
        else
            valer(7) = zr(jmax+6)
            symax = r8vide()
            call getvr8(' ', 'SY_MAX', scal=symax, nbret=n1)
            if (n1 .eq. 0) then
                call rcvale(mater, 'RCCM', 0, k8b, [rbid], &
                            1, 'SY_02   ', valres(1), icodre(1), 0)
                if (icodre(1) .eq. 0) then
                    symax = valres(1)
                else
                    call utmess('F', 'POSTRCCM_4')
                end if
            end if
            call rcmcrt(symax, sigmoypres, valer(8), valer(10))
        end if
!
        valer(9) = zr(jmax+7)
        do k = 1, 3
            valer(10+k) = zr(jmax-1+8+k)
        end do
        if (lefat) then
            call getvr8('ENVIRONNEMENT', 'FEN_INTEGRE', iocc=1, scal=fenint, nbret=n5)
            if (n5 .eq. 0 .or. abs(fenint) .lt. 1e-8) call utmess('F', 'POSTRCCM_54')
            fenglobal = 0.d0
            if (abs(zr(jmax+10)) .gt. 1e-8) fenglobal = zr(jmax+12)/zr(jmax+10)
            valer(16) = zr(jmax+10)
            if (fenglobal .gt. fenint) valer(16) = zr(jmax+12)/fenint
            valer(14) = fenglobal
            valer(15) = fenint
        else
            valer(14) = r8vide()
            valer(15) = r8vide()
        end if
!
!
        if (typtab .ne. 'SYSTUS') then
            call tbajli(nomres, npar1, nopar1, [ibid], valer, &
                        [c16b], valek, 0)
        end if
    end do
!
    if (typtab .eq. 'VALE_MAX') goto 999
!
!
!     ----------------------------------------------------------------
! --- AFFICHAGE DES GRANDEURS PAR SITUATION
!     ----------------------------------------------------------------
    call getfac('SITUATION', nb)
    call getfac('SEISME', ns)
    call jeveuo('&&RC3200.SITU_INFOI', 'L', jval)
    call jeveuo('&&RC3200.SITU_NOM', 'L', jnom)
!
    if (typtab .eq. 'SYSTUS') goto 444
!
!
    valek(1) = 'SITU'
!
    do i = 1, nb, 1
! --aller chercher le numéro de situation
        valei(1) = zi(jval+27*(i-1))
! aller chercher le numéro de groupe
        valei(2) = zi(jval+27*(i-1)+1)
! --aller chercher le nom de situation
        valek(4) = zk16(jnom+(i-1))
!
        do im = 1, 2
            call jeveuo('&&RC3200.SITU_RESU.'//lieu(im), 'L', jresu)
            valek(3) = lieu(im)
            valek(2) = 'SANS'
!
            do k = 1, 11
                valer(k) = zr(jresu+123*(i-1)-1+k)
            end do
!
            sigmoypres = zr(jresu+123*(i-1)+9)
            if (sigmoypres .eq. r8vide()) then
                valer(12) = r8vide()
                valer(14) = r8vide()
            else
                symax = r8vide()
                call getvr8(' ', 'SY_MAX', scal=symax, nbret=n1)
                if (n1 .eq. 0) then
                    call rcvale(mater, 'RCCM', 0, k8b, [rbid], &
                                1, 'SY_02   ', valres(1), icodre(1), 0)
                    if (icodre(1) .eq. 0) then
                        symax = valres(1)
                    else
                        call utmess('F', 'POSTRCCM_4')
                    end if
                end if
                call rcmcrt(symax, sigmoypres, valer(12), valer(13))
            end if
!
            valer(14) = zr(jresu+123*(i-1)-1+122)
            valer(15) = zr(jresu+123*(i-1)-1+123)
            valer(16) = zr(jresu+123*(i-1)-1+18)
            valer(17) = zr(jresu+123*(i-1)-1+17)
            valer(18) = zr(jresu+123*(i-1)-1+12)
            valer(19) = zr(jresu+123*(i-1)-1+13)
            valer(20) = zr(jresu+123*(i-1)-1+14)
            valer(21) = zr(jresu+123*(i-1)-1+15)
            if (sscyc .eq. 'OUI') then
                valer(22) = zr(jresu+123*(i-1)-1+121)
            else
                valer(22) = r8vide()
            end if
            valer(23) = zr(jresu+123*(i-1)-1+16)
!
            call tbajli(nomres, npar2, nopar2, valei, valer, &
                        [c16b], valek, 0)
!
            if (ns .eq. 0) goto 888
!
            call jeveuo('&&RC3200.SITUS_RESU.'//lieu(im), 'L', jresus)
            valek(2) = 'AVEC'
!
            do k = 1, 11
                valer(k) = zr(jresus+123*(i-1)-1+k)
            end do
!
            valer(14) = zr(jresus+123*(i-1)-1+122)
            valer(15) = zr(jresus+123*(i-1)-1+123)
            valer(16) = zr(jresus+123*(i-1)-1+18)
            valer(17) = zr(jresus+123*(i-1)-1+17)
            valer(18) = zr(jresus+123*(i-1)-1+12)
            valer(19) = zr(jresus+123*(i-1)-1+13)
            valer(20) = zr(jresus+123*(i-1)-1+14)
            valer(21) = zr(jresus+123*(i-1)-1+15)
            if (sscyc .eq. 'OUI') then
                valer(22) = zr(jresus+123*(i-1)-1+121)
            else
                valer(22) = r8vide()
            end if
            valer(23) = zr(jresus+123*(i-1)-1+16)
!
            call tbajli(nomres, npar2, nopar2, valei, valer, &
                        [c16b], valek, 0)
!
888         continue
!
        end do
    end do
!
    if (.not. lfat) goto 999
!
!
!     ----------------------------------------------------------------
! --- AFFICHAGE DES GRANDEURS PAR COMBINAISON
!     ----------------------------------------------------------------
! cette combinaison est-elle possible ?
    call jeveuo('&&RC3200.COMBI', 'L', jcombi)
    valek(1) = 'COMB'
!
    do iocc1 = 1, nb
        do iocc2 = 1, nb
            if (zi(jcombi+nb*(iocc1-1)+iocc2-1) .ne. 0 .and. iocc1 .ne. iocc2) then
! ----------- aller chercher les numéros de situation
                valei(1) = zi(jval+27*(iocc1-1))
                valei(3) = zi(jval+27*(iocc2-1))
! ----------- aller chercher les numéros de groupe
                valei(2) = zi(jval+27*(iocc1-1)+1)
                valei(4) = zi(jval+27*(iocc2-1)+1)
!
                do im = 1, 2
! ----------- aller chercher les noms de situation
                    valek(4) = zk16(jnom+(iocc1-1))
                    valek(5) = zk16(jnom+(iocc2-1))
                    call jeveuo('&&RC3200.COMB_RESU.'//lieu(im), 'L', jresucomb)
                    valek(3) = lieu(im)
                    valek(2) = 'SANS'
                    do k = 1, 6
                        valer(k) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+k)
                    end do
!
                    sigmoypres = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+7)
                    valer(7) = sigmoypres
                    if (sigmoypres .eq. r8vide()) then
                        valer(9) = r8vide()
                        valer(10) = r8vide()
                    else
                        symax = r8vide()
                        call getvr8(' ', 'SY_MAX', scal=symax, nbret=n1)
                        if (n1 .eq. 0) then
                            call rcvale(mater, 'RCCM', 0, k8b, [rbid], &
                                        1, 'SY_02   ', valres(1), icodre(1), 0)
                            if (icodre(1) .eq. 0) then
                                symax = valres(1)
                            else
                                call utmess('F', 'POSTRCCM_4')
                            end if
                        end if
                        call rcmcrt(symax, sigmoypres, valer(9), valer(10))
                    end if
!
!
                    valer(8) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+8)
                    valer(11) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+22)
                    valer(12) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+23)
!
                    do k = 1, 6
                        valer(12+k) = r8vide()
                    end do
!
                    if (sscyc .eq. 'OUI') then
                        valer(19) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+21)
                    else
                        valer(19) = r8vide()
                    end if
                    valer(20) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+17)
                    call tbajli(nomres, npar3, nopar3, valei, valer, &
                                [c16b], valek, 0)
! --------------- on crée la ligne avec le transitoires fictifs uniquement
                    valek(4) = 'FICTIF1'
                    valek(5) = 'FICTIF1'
!
                    do k = 1, 12
                        valer(k) = r8vide()
                    end do
!
                    valer(13) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+9)- &
                                zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+18)
                    valer(14) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+18)
                    valer(15) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+9)
                    valer(16) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+10)
                    valer(17) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+11)
                    valer(18) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+12)
                    valer(19) = r8vide()
                    valer(20) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+24)
!
                    call tbajli(nomres, npar3, nopar3, valei, valer, &
                                [c16b], valek, 0)
!
                    valek(4) = 'FICTIF2'
                    valek(5) = 'FICTIF2'
!
                    valer(13) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+13)- &
                                zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+19)
                    valer(14) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+19)
                    valer(15) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+13)
                    valer(16) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+14)
                    valer(17) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+15)
                    valer(18) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+16)
                    valer(19) = r8vide()
                    valer(20) = zr(jresucomb+25*nb*(iocc1-1)+25*(iocc2-1)-1+25)
!
                    call tbajli(nomres, npar3, nopar3, valei, valer, &
                                [c16b], valek, 0)
!
                    if (ns .eq. 0) goto 777
!
                    valek(2) = 'AVEC'
                    valek(4) = zk16(jnom+(iocc1-1))
                    valek(5) = zk16(jnom+(iocc2-1))
                    call jeveuo('&&RC3200.COMBS_RESU.'//lieu(im), 'L', jresucombs)
                    valek(3) = lieu(im)
                    do k = 1, 6
                        valer(k) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+k)
                    end do
!
                    sigmoypres = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+7)
                    valer(7) = sigmoypres
                    if (sigmoypres .eq. r8vide()) then
                        valer(9) = r8vide()
                        valer(10) = r8vide()
                    else
                        symax = r8vide()
                        call getvr8(' ', 'SY_MAX', scal=symax, nbret=n1)
                        if (n1 .eq. 0) then
                            call rcvale(mater, 'RCCM', 0, k8b, [rbid], &
                                        1, 'SY_02   ', valres(1), icodre(1), 0)
                            if (icodre(1) .eq. 0) then
                                symax = valres(1)
                            else
                                call utmess('F', 'POSTRCCM_4')
                            end if
                        end if
                        call rcmcrt(symax, sigmoypres, valer(9), valer(10))
                    end if
!
                    valer(8) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+8)
                    valer(11) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+22)
                    valer(12) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+23)
                    do k = 1, 6
                        valer(12+k) = r8vide()
                    end do
                    if (sscyc .eq. 'OUI') then
                        valer(19) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+21)
                    else
                        valer(19) = r8vide()
                    end if
                    valer(20) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+17)
!
                    call tbajli(nomres, npar3, nopar3, valei, valer, &
                                [c16b], valek, 0)
! --------------- on crée la ligne avec le transitoires fictifs uniquement
                    valek(4) = 'FICTIF1'
                    valek(5) = 'FICTIF1'
!
                    do k = 1, 12
                        valer(k) = r8vide()
                    end do
                    valer(13) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+9)- &
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+18)
                    valer(14) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+18)
                    valer(15) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+9)
                    valer(16) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+10)
                    valer(17) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+11)
                    valer(18) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+12)
                    valer(19) = r8vide()
                    valer(20) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+24)
!
                    call tbajli(nomres, npar3, nopar3, valei, valer, &
                                [c16b], valek, 0)
!
                    valek(4) = 'FICTIF2'
                    valek(5) = 'FICTIF2'
!
                    valer(13) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+13)- &
                                zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+19)
                    valer(14) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+19)
                    valer(15) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+13)
                    valer(16) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+14)
                    valer(17) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+15)
                    valer(18) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+16)
                    valer(19) = r8vide()
                    valer(20) = zr(jresucombs+25*nb*(iocc1-1)+25*(iocc2-1)-1+25)
!
                    call tbajli(nomres, npar3, nopar3, valei, valer, &
                                [c16b], valek, 0)
777                 continue
!
                end do
            end if
        end do
    end do
444 continue
!
!     ----------------------------------------------------------------
! --- AFFICHAGE DES GRANDEURS QUI INTERVIENNENT DANS FU_TOTAL
! --- SN,'INST_SN_1','INST_SN_2', 'SN*', 'INST_SN*_1', 'INST_SN*_2',
! --- 'KE_MECA', 'KE_THER','SP(MECA)', 'SP(THER)', 'SP',
! --- 'INST_SALT_1', 'INST_SALT_2','SALT', 'FU_UNIT_SS_CYCL', 'FU_UNIT', 'NOCC_PRIS', 'FU_PARTIEL'
! --- , 'FEN', 'FEN_ELAS','FUEN_PARTIEL'
!     ----------------------------------------------------------------
!
    valek(1) = 'FACT'
!
    do im = 1, 2
        valek(3) = lieu(im)
        call jeveuo('&&RC3200.FACT.'//lieu(im), 'L', jfact)
        if (lefat) call jeveuo('&&RC3200.FACTENV.'//lieu(im), 'L', jfactenv)
!
        k = 0
!
555     continue
!
        num1 = zi(jfact+6*k)
        num2 = zi(jfact+6*k+1)
!
        if (num1 .eq. 0) goto 666
        noccpris = zi(jfact+6*k+4)
!
        if (zi(jfact+6*k+5) .eq. 2) then
            valek(2) = 'AVEC'
            call jeveuo('&&RC3200.SITUS_RESU.'//lieu(im), 'L', ind1)
            call jeveuo('&&RC3200.COMBS_RESU.'//lieu(im), 'L', ind2)
            call jeveuo('&&RC3200.MAX_RESU.'//lieu(im), 'L', jmax)
            fuseism = zr(jmax+11)
        end if
        if (zi(jfact+6*k+5) .eq. 1) then
            valek(2) = 'SANS'
            call jeveuo('&&RC3200.SITU_RESU.'//lieu(im), 'L', ind1)
            call jeveuo('&&RC3200.COMB_RESU.'//lieu(im), 'L', ind2)
            fuseism = 0.d0
        end if
!
        valek(4) = zk16(jnom+(num1-1))
        valek(5) = zk16(jnom+(num2-1))
!
!---- numéro de situation
        valei(1) = zi(jval+27*(num1-1))
!---- nombre d'occurences restantes de la situation
        valei(2) = zi(jfact+6*k+2)
!---- groupe de la situation
        valei(3) = zi(jval+27*(num1-1)+1)
!
!---- numéro de situation
        valei(4) = zi(jval+27*(num2-1))
!---- nombre d'occurences restantes de la situation
        valei(5) = zi(jfact+6*k+3)
!---- groupe de la situation
        valei(6) = zi(jval+27*(num2-1)+1)
!
        valei(7) = noccpris
!
!---- une situation seule a le plus grand fu unitaire
        if (num1 .eq. num2) then
!
            valer(1) = zr(ind1+123*(num1-1)+3)
            valer(2) = zr(ind1+123*(num1-1)+4)
            valer(3) = zr(ind1+123*(num1-1)+5)
            valer(4) = zr(ind1+123*(num1-1)+6)
            valer(5) = zr(ind1+123*(num1-1)+7)
            valer(6) = zr(ind1+123*(num1-1)+8)
            valer(7) = zr(ind1+123*(num1-1)+121)
            valer(8) = zr(ind1+123*(num1-1)+122)
            valer(9) = zr(ind1+123*(num1-1)+17)
            valer(10) = zr(ind1+123*(num1-1)+16)
            valer(11) = zr(ind1+123*(num1-1)+11)
!
            valer(12) = zr(ind1+123*(num1-1)+12)
            valer(13) = zr(ind1+123*(num1-1)+13)
            valer(14) = zr(ind1+123*(num1-1)+14)
            if (sscyc .eq. 'OUI') then
                valer(15) = zr(ind1+123*(num1-1)+120)
            else
                valer(15) = r8vide()
            end if
!
            valer(16) = zr(ind1+123*(num1-1)+15)+fuseism
            valer(17) = valer(16)*noccpris
!
            if (lefat) then
                valer(18) = zr(jfactenv+3*k)
                valer(19) = zr(jfactenv+3*k+1)
                valer(20) = zr(jfactenv+3*k+2)
            else
                valer(18) = r8vide()
                valer(19) = r8vide()
                valer(20) = r8vide()
            end if
!
            call tbajli(nomres, npar4, nopar4, valei, valer, &
                        [c16b], valek, 0)
!
!---- une combinaison de situations a le plus grand fu unitaire
        else
!
            valer(1) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+1)
            valer(2) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+2)
            valer(3) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+3)
            valer(4) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+4)
            valer(5) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+5)
            valer(6) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+6)
            valer(7) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+22)
            valer(8) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+23)
            valer(9) = r8vide()
            valer(10) = r8vide()
            valer(11) = r8vide()
            valer(12) = r8vide()
            valer(13) = r8vide()
            valer(14) = r8vide()
            if (sscyc .eq. 'OUI') then
                valer(15) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+21)
            else
                valer(15) = r8vide()
            end if
            valer(16) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+17)+fuseism
            valer(17) = valer(16)*noccpris
!
            if (lefat) then
                valer(18) = zr(jfactenv+3*k)
                valer(19) = zr(jfactenv+3*k+1)
                valer(20) = zr(jfactenv+3*k+2)
            else
                valer(18) = r8vide()
                valer(19) = r8vide()
                valer(20) = r8vide()
            end if
!
!
            call tbajli(nomres, npar4, nopar4, valei, valer, &
                        [c16b], valek, 0)
!
            valek(4) = 'FICTIF1'
            valek(5) = 'FICTIF1'
            do i = 1, 8
                valer(i) = r8vide()
            end do
!
            valer(9) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+9)- &
                       zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+18)
            valer(10) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+18)
            valer(11) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+9)
            valer(12) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+10)
            valer(13) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+11)
            valer(14) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+12)
            valer(15) = r8vide()
            valer(16) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+24)
            valer(17) = r8vide()
            valer(18) = r8vide()
            valer(19) = r8vide()
            valer(20) = r8vide()
!
            call tbajli(nomres, npar4, nopar4, valei, valer, &
                        [c16b], valek, 0)
!
            valek(4) = 'FICTIF2'
            valek(5) = 'FICTIF2'
!
!
            valer(9) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+13)- &
                       zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+19)
            valer(10) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+19)
            valer(11) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+13)
            valer(12) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+14)
            valer(13) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+15)
            valer(14) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+16)
            valer(16) = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+25)
!
            call tbajli(nomres, npar4, nopar4, valei, valer, &
                        [c16b], valek, 0)
!
        end if
!
        k = k+1
        goto 555
!
666     continue
!
    end do
!
    if (typtab .eq. 'SYSTUS') then
!
        do i = 1, 5
            valek(i) = ' '
        end do
!
        do i = 1, 7
            valei(i) = 0
        end do
!
        valek(1) = 'TOTAL'
        do im = 1, 2
            valek(3) = lieu(im)
!
            call jeveuo('&&RC3200.MAX_RESU.'//lieu(im), 'L', jmax)
!
            do i = 1, 20
                valer(i) = r8vide()
            end do
!
            valer(17) = zr(jmax-1+11)
            valer(20) = zr(jmax+12)
            if (lefat) then
                call getvr8('ENVIRONNEMENT', 'FEN_INTEGRE', iocc=1, scal=fenint, nbret=n5)
                if (n5 .eq. 0 .or. abs(fenint) .lt. 1e-8) call utmess('F', 'POSTRCCM_54')
                fenglobal = 0.d0
                if (abs(zr(jmax+10)) .gt. 1e-8) fenglobal = zr(jmax+12)/zr(jmax+10)
                if (fenglobal .gt. fenint) valer(17) = zr(jmax+12)/fenint
            end if
!
            call tbajli(nomres, npar4, nopar4, valei, valer, &
                        [c16b], valek, 0)
        end do
    end if
!
999 continue
!
    call jedema()
!
end subroutine
