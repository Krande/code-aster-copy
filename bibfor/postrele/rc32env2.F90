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
subroutine rc32env2(iocc1, iocc2, ke_pond, lieu, fen)
    implicit none
!
    integer(kind=8) :: iocc1, iocc2
    real(kind=8) :: ke_pond, fen
    character(len=4) :: lieu
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE B3200, ZE200
!     CALCUL DU FACTEUR D'ENVIRONNEMENT D'UNE COMBINAISON DE TRANSITOIRES :
!        IN  : iocc1    :  NUMERO DU PREMIER TRANSITOIRE (SITUATION iocc1)
!        IN  : iocc2    :  NUMERO DU DEUXIEME TRANSITOIRE (SITUATION iocc2)
!        IN  : KE_POND   :  KE PONDERE DE LA COMBINAISON DE CES 2 TRANSITOIRES
!        IN  : LIEU :  ORIG OU EXTR
!        OUT : FEN  :  FACTEUR D'ENVIRONNEMENT DE LA COMBINAISON
!                      SI KE_POND = 1.0, FEN ELASTIQUE
!
!     ------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lcqeqv.h"
#include "asterfort/rcjaco.h"
#include "asterfort/rcveri.h"
#include "asterfort/tbexv1.h"
#include "asterfort/tbliva.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=16) :: motclf, motclf2, motclf3, motclf4, valek(2)
    character(len=16) :: motclf5, val
    integer(kind=8) :: nbther, nbpres, nbmeca, situ(2), ncmp
    integer(kind=8) :: n1, n2, n3, nume1, nume2, nume3, ither, ipres, imeca
    integer(kind=8) :: numether, numepres, numemeca, n4, nbinst, jinst, nbabsc, jabsc
    integer(kind=8) :: i, j, m, ibid, iret, ndim, jcont, nbinst2, jdepsi, itrouve
    integer(kind=8) :: jj, jfen, n5(14), jtemp, n6, n7, jtempy, nbtempy, ii
    parameter(ncmp=6)
    character(len=8) :: nocmp(ncmp), crit(2), table1, table2, table3
    character(len=8) :: tableok, k8b, table4, table5
    real(kind=8) :: prec(2), vale(2), diff(6), nul(6), equi(3)
    real(kind=8) :: tresca, contraintesth, contraintespr, contraintesmec
    real(kind=8) :: prinmax, e, nume(2), deno(2), epset, set, oet, tet
    real(kind=8) :: a, b, c, epsisup, epsiinf, temp, tmoy, tsup, tinf
    real(kind=8) :: valtinf, valtsup, valtmoynum, valtmoyden, epsiet
    real(kind=8) :: critepsi, tempmin, tempmax, tempii, emin, emax
    character(len=24) :: instan, abscur, instany, valk(4)
    complex(kind=8) :: cbid
!
!
! DEB ------------------------------------------------------------------
    call jemarq()
!
!  ------ POUR LES DEUX SITUATION iocc1 ET iocc2
!  ------ ON STOCKE LE TRANSITOIRE ASSOCIE
!
    motclf = 'RESU_THER'
    motclf2 = 'RESU_PRES'
    motclf3 = 'RESU_MECA'
    motclf4 = 'SITUATION'
    motclf5 = 'ENVIRONNEMENT'
!
    call getfac(motclf, nbther)
    call getfac(motclf2, nbpres)
    call getfac(motclf3, nbmeca)
!
    nocmp(1) = 'SIXX'
    nocmp(2) = 'SIXY'
    nocmp(3) = 'SIXZ'
    nocmp(4) = 'SIYY'
    nocmp(5) = 'SIYZ'
    nocmp(6) = 'SIZZ'
!
    valek(1) = 'INST            '
    valek(2) = 'ABSC_CURV       '
    val = 'TEMP       '
!
    prec(1) = 1.0d-06
    prec(2) = 1.0d-06
    crit(1) = 'RELATIF'
    crit(2) = 'RELATIF'
!
    situ(1) = iocc1
    situ(2) = iocc2
!
    nul(1) = 0.d0
    nul(2) = 0.d0
    nul(3) = 0.d0
    nul(4) = 0.d0
    nul(5) = 0.d0
    nul(6) = 0.d0
!
    call getvr8(motclf5, 'S_ETOILE', iocc=1, scal=set, nbret=n5(1))
    call getvr8(motclf5, 'SEUIL_EPSI_SUP', iocc=1, scal=epsisup, nbret=n5(2))
    call getvr8(motclf5, 'SEUIL_EPSI_INF', iocc=1, scal=epsiinf, nbret=n5(3))
    call getvr8(motclf5, 'A_ENV', iocc=1, scal=a, nbret=n5(4))
    call getvr8(motclf5, 'B_ENV', iocc=1, scal=b, nbret=n5(5))
    call getvr8(motclf5, 'C_ENV', iocc=1, scal=c, nbret=n5(6))
    call getvr8(motclf5, 'SEUIL_T_SUP', iocc=1, scal=tsup, nbret=n5(7))
    call getvr8(motclf5, 'SEUIL_T_INF', iocc=1, scal=tinf, nbret=n5(8))
    call getvr8(motclf5, 'VALE_T_INF', iocc=1, scal=valtinf, nbret=n5(9))
    call getvr8(motclf5, 'VALE_T_SUP', iocc=1, scal=valtsup, nbret=n5(10))
    call getvr8(motclf5, 'VALE_T_MOY_NUM', iocc=1, scal=valtmoynum, nbret=n5(11))
    call getvr8(motclf5, 'VALE_T_MOY_DEN', iocc=1, scal=valtmoyden, nbret=n5(12))
    call getvr8(motclf5, 'CRIT_EPSI', iocc=1, scal=critepsi, nbret=n5(13))
    critepsi = critepsi/100.d0
!
    call getvid(motclf5, 'TABL_YOUNG', iocc=1, scal=table5, nbret=n5(14))
!
    do i = 1, 14
        if (n5(i) .eq. 0) call utmess('F', 'POSTRCCM_54')
    end do
!
    instany = '&&RC32.INSTANTYOUNG'
    call tbexv1(table5, 'TEMP', instany, 'V', nbtempy, &
                k8b)
    call jeveuo(instany, 'L', jtempy)
!
    do i = 1, 2
!
        contraintesth = 0.d0
        contraintespr = 0.d0
        contraintesmec = 0.d0
!
        itrouve = situ(i)
!
        call getvr8(motclf4, 'O_ETOILE', iocc=itrouve, scal=oet, nbret=n6)
        call getvid(motclf4, 'TABL_TEMP', iocc=itrouve, scal=table4, nbret=n7)
!
        call getvis(motclf4, 'NUME_RESU_THER', iocc=itrouve, scal=nume1, nbret=n1)
        call getvis(motclf4, 'NUME_RESU_PRES', iocc=itrouve, scal=nume2, nbret=n2)
        call getvis(motclf4, 'NUME_RESU_MECA', iocc=itrouve, scal=nume3, nbret=n3)
!
        if (n1 .ne. 0) then
            do ither = 1, nbther, 1
                call getvis(motclf, 'NUME_RESU_THER', iocc=ither, scal=numether, nbret=n4)
                if (numether .eq. nume1) then
                    call getvid(motclf, 'TABL_RESU_THER', iocc=ither, scal=table1, nbret=n4)
                end if
            end do
        end if
!
        if (n2 .ne. 0) then
            do ipres = 1, nbpres, 1
                call getvis(motclf2, 'NUME_RESU_PRES', iocc=ipres, scal=numepres, nbret=n4)
                if (numepres .eq. nume2) then
                    call getvid(motclf2, 'TABL_RESU_PRES', iocc=ipres, scal=table2, nbret=n4)
                end if
            end do
        end if
!
        if (n3 .ne. 0) then
            do imeca = 1, nbmeca, 1
                call getvis(motclf3, 'NUME_RESU_MECA', iocc=imeca, scal=numemeca, nbret=n4)
                if (numemeca .eq. nume3) then
                    call getvid(motclf3, 'TABL_RESU_MECA', iocc=imeca, scal=table3, nbret=n4)
                end if
            end do
        end if
!
! ------ RECUPERATION DES INSTANTS ET DES ABSCISSES
!
! --------- grace a la table thermique ou de pression ou mecanique
!
        if (n1 .ne. 0) then
            tableok = table1
        else if (n2 .ne. 0) then
            tableok = table2
        else
            tableok = table3
        end if
!
! --------- on verifie l'ordre des noeuds de la table
!
        call rcveri(tableok)
!
! --------- on recupere les instants de la table
!
        instan = '&&RC32.INSTANT'
        call tbexv1(tableok, valek(1), instan, 'V', nbinst, &
                    k8b)
        call jeveuo(instan, 'L', jinst)
!
! --------- on recupere les abscisses curvilignes de la table
!
        abscur = '&&RC32.ABSC_CURV'
        call tbexv1(tableok, valek(2), abscur, 'V', nbabsc, &
                    k8b)
        call jeveuo(abscur, 'L', jabsc)
!
        if (lieu .eq. 'ORIG') then
            vale(2) = zr(jabsc)
        else
            vale(2) = zr(jabsc+nbabsc-1)
        end if
!
! ------ STOCKAGE DES CONTRAINTES POUR LE CALCUL DU FEN
!
        ndim = ncmp*nbinst
        nbinst2 = nbinst-1
        call wkvect('&&RC32.CONT', 'V V R', ndim, jcont)
        call wkvect('&&RC32.DEPSI', 'V V R', nbinst2, jdepsi)
        call wkvect('&&RC32.FENJ', 'V V R', nbinst2, jfen)
        call wkvect('&&RC32.TEMP', 'V V R', nbinst2, jtemp)
!
        do j = 1, nbinst
            vale(1) = zr(jinst+j-1)
!
! --------- on recupere la temperature a l 'instant i
!
            call tbliva(table4, 2, valek, [ibid], vale, &
                        [cbid], k8b, crit, prec, 'TEMP', &
                        k8b, ibid, temp, cbid, k8b, &
                        iret)
            if (iret .ne. 0) then
                valk(1) = table4
                valk(2) = 'TEMP'
                valk(3) = valek(1)
                valk(4) = valek(2)
                call utmess('F', 'POSTRCCM_2', nk=4, valk=valk, nr=2, &
                            valr=vale)
            end if
            zr(jtemp+j-1) = temp
!
            do m = 1, ncmp
!
! --------- on recupere les contraintes(t) a l 'instant i
!
                if (n1 .ne. 0) then
                    call tbliva(table1, 2, valek, [ibid], vale, &
                                [cbid], k8b, crit, prec, nocmp(m), &
                                k8b, ibid, contraintesth, cbid, k8b, &
                                iret)
                end if
!
                if (n2 .ne. 0) then
                    call tbliva(table2, 2, valek, [ibid], vale, &
                                [cbid], k8b, crit, prec, nocmp(m), &
                                k8b, ibid, contraintespr, cbid, k8b, &
                                iret)
                end if
!
                if (n3 .ne. 0) then
                    call tbliva(table3, 2, valek, [ibid], vale, &
                                [cbid], k8b, crit, prec, nocmp(m), &
                                k8b, ibid, contraintesmec, cbid, k8b, &
                                iret)
                end if
                zr(jcont+(j-1)*ncmp+m-1) = contraintesth+contraintespr+contraintesmec
                if (j .ge. 2) then
                    diff(m) = zr(jcont+(j-1)*ncmp+m-1)-zr(jcont+(j-2)*ncmp+m-1)
                    if (m .eq. ncmp) then
                        tmoy = (zr(jtemp+j-1)+zr(jtemp+j-2))/2
                        if (lcqeqv(diff, nul) .eq. 'OUI') then
                            tresca = 0.d0
                        else
                            call rcjaco(diff, equi)
! --------- on recupere cherche la contrainte principale maxi
                            if (abs(equi(1)) .gt. abs(equi(2))) then
                                prinmax = equi(1)
                            else
                                prinmax = equi(2)
                            end if
                            if (abs(prinmax) .lt. abs(equi(3))) then
                                prinmax = equi(3)
                            end if
! --------- si cette cont. princ. max est <0, la vitesse de déformation est nulle
                            if (prinmax .lt. 0) then
                                tresca = 0.d0
! --------- si cette cont. princ. max est >0,
! --------- la vitesse de déformation vaut ke_pond*tresca/(E*(ti+1-ti))
                            else
                                tresca = max( &
                                         abs(equi(1)-equi(2)), abs(equi(1)-equi(3)), &
                                         abs(equi(2)-equi(3)) &
                                         )
                            end if
                        end if
! --------- calcul de e
                        tempmin = zr(jtempy)
                        tempmax = zr(jtempy)
                        ii = 0
120                     continue
                        ii = ii+1
                        tempii = zr(jtempy+ii)
                        if (tmoy .ge. tempii) then
                            tempmin = tempii
                            goto 120
                        end if
                        if (abs(tmoy-tempmin) .lt. 1e-6) then
                            tempmax = tmoy
                        else
                            tempmax = tempii
                        end if
                        call tbliva(table5, 1, val, [ibid], [tempmin], &
                                    [cbid], k8b, crit, prec, 'YOUNG', &
                                    k8b, ibid, emin, cbid, k8b, &
                                    iret)
                        if (iret .eq. 1) then
                            valk(1) = table5
                            valk(2) = 'YOUNG'
                            call utmess('F', 'POSTRCCM_1', nk=2, valk=valk)
                        end if
                        call tbliva(table5, 1, val, [ibid], [tempmax], &
                                    [cbid], k8b, crit, prec, 'YOUNG', &
                                    k8b, ibid, emax, cbid, k8b, &
                                    iret)
                        if (abs(tempmin-tempmax) .lt. 1e-6) then
                            e = emin
                        else
                            e = ((emin-emax)/(tempmin-tempmax))*(tmoy-tempmin)+emin
                        end if
! --------- calcul de delta_epsilon
                        zr(jdepsi+j-2) = ke_pond*tresca/e
! --------- calcul de epsilon_point puis epsilon_point*
                        epsiet = zr(jdepsi+j-2)/(zr(jinst+j-1)-zr(jinst+j-2))
                        if (epsiet*100 .lt. epsiinf) then
                            epset = log10(epsiinf/epsisup)/log10(exp(1.d0))
                        else if (epsiet*100 .gt. epsisup) then
                            epset = 0.d0
                        else
                            epset = log10(100*epsiet/epsisup)/log10(exp(1.d0))
                        end if
! --------- calcul de t*
                        if (tmoy .lt. tinf) then
                            tet = valtinf
                        else if (tmoy .gt. tsup) then
                            tet = valtsup
                        else
                            tet = (tmoy-valtmoynum)/valtmoyden
                        end if
!
                        zr(jfen+j-2) = exp((a+b*epset)*set*oet*tet+c)
                    end if
                end if
            end do
        end do
!
        nume(i) = 0.d0
        deno(i) = 0.d0
        do jj = 1, nbinst-1
            nume(i) = nume(i)+(zr(jdepsi+jj-1))*(zr(jfen+jj-1))
            deno(i) = deno(i)+zr(jdepsi+jj-1)
        end do
!
        call jedetr(instan)
        call jedetr(abscur)
        call jedetr('&&RC32.CONT')
        call jedetr('&&RC32.DEPSI')
        call jedetr('&&RC32.FENJ')
        call jedetr('&&RC32.TEMP')
!
    end do
!
    call jedetr(instany)
!
    if (deno(1)+deno(2) .le. critepsi) then
        fen = 1
    else
        fen = (nume(1)+nume(2))/(deno(1)+deno(2))
    end if
!
    call jedema()
end subroutine
