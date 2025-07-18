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

subroutine rvmima(nomres, iocc)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/prexel.h"
#include "asterfort/prexno.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnopa.h"
#include "asterfort/rsutnu.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbexip.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: iocc
    character(len=*) :: nomres
!
!     COMMANDE : POST_RELEVE, OPERATION='EXTREMA'
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: nbpano, nbpael
    parameter(nbpano=6, nbpael=7)
    character(len=16) :: nopano(nbpano), nopael(nbpael), nopara(200)
!
    integer(kind=8) :: ibid, n1, np, nc, iret
    integer(kind=8) :: jordr, i100, nbordr, iord, vali(20), nbc, nbpar
    integer(kind=8) :: ispmax, ispmin, isamax, isamin
    integer(kind=8) :: iadr, iac, ii, ik, ir, jaces, nbacc
    real(kind=8) :: prec, valr(200), valmax, valmin, vaamax, vaamin
    complex(kind=8) :: c16b
    character(len=3) :: typpar
    character(len=8) :: crit, resu, mamax, nomax, mamin, nomin, tych, ctype
    character(len=8) :: maamax, noamax, maamin, noamin
    character(len=8) :: cpmax, cpmin, cpamax, cpamin
    character(len=16) :: nomcha, intitu
    character(len=19) :: knum, champ
    character(len=24) :: nomjv
    character(len=80) :: valk(200)
    aster_logical :: exist
!
    data nopano/'INTITULE', 'CHAM_GD',&
     &              'EXTREMA', 'NOEUD', 'CMP', 'VALE'/
    data nopael/'INTITULE', 'CHAM_GD',&
     &              'EXTREMA', 'MAILLE', 'NOEUD', 'CMP', 'VALE'/
! ---------------------------------------------------------------------
!
    call jemarq()
    c16b = (0.d0, 0.d0)
    knum = '&&RVMIMA.NUME_ORDRE'
    nbc = 0
!
    call getvtx('ACTION', 'INTITULE', iocc=iocc, scal=intitu, nbret=n1)
    nopara(1) = 'INTITULE'
    valk(1) = intitu
!
! ----- TRAITEMENT DU CHAMP_GD  -----
!
    call getvid('ACTION', 'CHAM_GD', iocc=iocc, scal=champ, nbret=n1)
    if (n1 .ne. 0) then
        valk(2) = champ
        call dismoi('TYPE_CHAMP', champ, 'CHAMP', repk=tych)
        if (tych(1:4) .eq. 'NOEU') then
            call prexno(champ, iocc, nomax, cpmax, valmax, &
                        nomin, cpmin, valmin, noamax, cpamax, &
                        vaamax, noamin, cpamin, vaamin)
            valr(1) = valmax
            valk(3) = 'MAX'
            valk(4) = nomax
            valk(5) = cpmax
            call tbajli(nomres, nbpano, nopano, vali, valr, &
                        [c16b], valk, 0)
            valr(1) = valmin
            valk(3) = 'MIN'
            valk(4) = nomin
            valk(5) = cpmin
            call tbajli(nomres, nbpano, nopano, vali, valr, &
                        [c16b], valk, 0)
            valr(1) = vaamax
            valk(3) = 'MAXI_ABS'
            valk(4) = noamax
            valk(5) = cpamax
            call tbajli(nomres, nbpano, nopano, vali, valr, &
                        [c16b], valk, 0)
            valr(1) = vaamin
            valk(3) = 'MINI_ABS'
            valk(4) = noamin
            valk(5) = cpamin
            call tbajli(nomres, nbpano, nopano, vali, valr, &
                        [c16b], valk, 0)
        else if (tych(1:4) .eq. 'ELNO') then
            call prexel(champ, iocc, mamax, nomax, ispmax, &
                        cpmax, valmax, mamin, nomin, ispmin, &
                        cpmin, valmin, maamax, noamax, isamax, &
                        cpamax, vaamax, maamin, noamin, isamin, &
                        cpamin, vaamin)
!
            valr(1) = valmax
            valk(3) = 'MAX'
            valk(4) = mamax
            valk(5) = nomax
            valk(6) = cpmax
            call tbajli(nomres, nbpael, nopael, vali, valr, &
                        [c16b], valk, 0)
            valr(1) = valmin
            valk(3) = 'MIN'
            valk(4) = mamin
            valk(5) = nomin
            valk(6) = cpmin
            call tbajli(nomres, nbpael, nopael, vali, valr, &
                        [c16b], valk, 0)
            valr(1) = vaamax
            valk(3) = 'MAXI_ABS'
            valk(4) = maamax
            valk(5) = noamax
            valk(6) = cpamax
            call tbajli(nomres, nbpael, nopael, vali, valr, &
                        [c16b], valk, 0)
            valr(1) = vaamin
            valk(3) = 'MINI_ABS'
            valk(4) = maamin
            valk(5) = noamin
            valk(6) = cpamin
            call tbajli(nomres, nbpael, nopael, vali, valr, &
                        [c16b], valk, 0)
!
        else
            call utmess('F', 'ALGORITH10_56', sk=tych)
        end if
        goto 999
    end if
!
! ----- TRAITEMENT DU RESULTAT  -----
!
    call getvid('ACTION', 'RESULTAT', iocc=iocc, scal=resu, nbret=n1)
    nopara(2) = 'RESU'
    valk(2) = resu
!
    call getvr8('ACTION', 'PRECISION', iocc=iocc, scal=prec, nbret=np)
    call getvtx('ACTION', 'CRITERE', iocc=iocc, scal=crit, nbret=nc)
    call rsutnu(resu, 'ACTION', iocc, knum, nbordr, &
                prec, crit, iret)
    if (iret .eq. 10) then
        call utmess('F', 'CALCULEL4_8', sk=resu)
    end if
    if (iret .ne. 0) then
        call utmess('F', 'ALGORITH3_41')
    end if
    call jeveuo(knum, 'L', jordr)
!
    call getvtx('ACTION', 'NOM_CHAM', iocc=iocc, scal=nomcha, nbret=nbc)
    nopara(3) = 'NOM_CHAM'
    valk(3) = nomcha
!
    do i100 = 1, nbordr
        iord = zi(jordr+i100-1)
!
        ik = 3
        ii = 0
        ir = 0
        nbpar = 3
!
        nbpar = nbpar+1
        nopara(nbpar) = 'NUME_ORDRE'
        ii = ii+1
        vali(ii) = iord
        nomjv = '&&RVMIMA.NOMS_ACCES'
        call rsnopa(resu, 0, nomjv, nbacc, ibid)
        if (nbacc .ne. 0) then
            call jeveuo(nomjv, 'L', jaces)
            do iac = 1, nbacc
                call rsadpa(resu, 'L', 1, zk16(jaces-1+iac), iord, &
                            1, sjv=iadr, styp=ctype, istop=0)
                call tbexip(nomres, zk16(jaces-1+iac), exist, typpar)
                if (.not. exist) then
                    call tbajpa(nomres, 1, zk16(jaces-1+iac), ctype)
                end if
                nbpar = nbpar+1
                nopara(nbpar) = zk16(jaces-1+iac)
                if (ctype(1:1) .eq. 'I') then
                    ii = ii+1
                    vali(ii) = zi(iadr)
                else if (ctype(1:1) .eq. 'R') then
                    ir = ir+1
                    valr(ir) = zr(iadr)
                else if (ctype(1:3) .eq. 'K80') then
                    ik = ik+1
                    valk(ik) = zk80(iadr)
                else if (ctype(1:3) .eq. 'K32') then
                    ik = ik+1
                    valk(ik) = zk32(iadr)
                else if (ctype(1:3) .eq. 'K24') then
                    ik = ik+1
                    valk(ik) = zk24(iadr)
                else if (ctype(1:3) .eq. 'K16') then
                    ik = ik+1
                    valk(ik) = zk16(iadr)
                else if (ctype(1:2) .eq. 'K8') then
                    ik = ik+1
                    valk(ik) = zk8(iadr)
                end if
            end do
            call jedetr(nomjv)
        end if
!
        nbpar = nbpar+1
        nopara(nbpar) = 'EXTREMA'
!
        call rsexch(' ', resu, nomcha, iord, champ, &
                    iret)
        if (iret .ne. 0) goto 100
        call dismoi('TYPE_CHAMP', champ, 'CHAMP', repk=tych)
!
        if (tych(1:4) .eq. 'NOEU') then
            nbpar = nbpar+1
            nopara(nbpar) = 'NOEUD'
            nbpar = nbpar+1
            nopara(nbpar) = 'CMP'
            nbpar = nbpar+1
            nopara(nbpar) = 'VALE'
!
            call prexno(champ, iocc, nomax, cpmax, valmax, &
                        nomin, cpmin, valmin, noamax, cpamax, &
                        vaamax, noamin, cpamin, vaamin)
            valr(ir+1) = valmax
            valk(ik+1) = 'MAX'
            valk(ik+2) = nomax
            valk(ik+3) = cpmax
            call tbajli(nomres, nbpar, nopara, vali, valr, &
                        [c16b], valk, 0)
            valr(ir+1) = valmin
            valk(ik+1) = 'MIN'
            valk(ik+2) = nomin
            valk(ik+3) = cpmin
            call tbajli(nomres, nbpar, nopara, vali, valr, &
                        [c16b], valk, 0)
            valr(ir+1) = vaamax
            valk(ik+1) = 'MAXI_ABS'
            valk(ik+2) = noamax
            valk(ik+3) = cpamax
            call tbajli(nomres, nbpar, nopara, vali, valr, &
                        [c16b], valk, 0)
            valr(ir+1) = vaamin
            valk(ik+1) = 'MINI_ABS'
            valk(ik+2) = noamin
            valk(ik+3) = cpamin
            call tbajli(nomres, nbpar, nopara, vali, valr, &
                        [c16b], valk, 0)
!
        else if (tych(1:4) .eq. 'ELNO') then
            nbpar = nbpar+1
            nopara(nbpar) = 'MAILLE'
            nbpar = nbpar+1
            nopara(nbpar) = 'NOEUD'
            nbpar = nbpar+1
            nopara(nbpar) = 'CMP'
            nbpar = nbpar+1
            nopara(nbpar) = 'VALE'
!
            call prexel(champ, iocc, mamax, nomax, ispmax, &
                        cpmax, valmax, mamin, nomin, ispmin, &
                        cpmin, valmin, maamax, noamax, isamax, &
                        cpamax, vaamax, maamin, noamin, isamin, &
                        cpamin, vaamin)
!
            valr(ir+1) = valmax
            valk(ik+1) = 'MAX'
            valk(ik+2) = mamax
            valk(ik+3) = nomax
            valk(ik+4) = cpmax
            call tbajli(nomres, nbpar, nopara, vali, valr, &
                        [c16b], valk, 0)
            valr(ir+1) = valmin
            valk(ik+1) = 'MIN'
            valk(ik+2) = mamin
            valk(ik+3) = nomin
            valk(ik+4) = cpmin
            call tbajli(nomres, nbpar, nopara, vali, valr, &
                        [c16b], valk, 0)
            valr(ir+1) = vaamax
            valk(ik+1) = 'MAXI_ABS'
            valk(ik+2) = maamax
            valk(ik+3) = noamax
            valk(ik+4) = cpamax
            call tbajli(nomres, nbpar, nopara, vali, valr, &
                        [c16b], valk, 0)
            valr(ir+1) = vaamin
            valk(ik+1) = 'MINI_ABS'
            valk(ik+2) = maamin
            valk(ik+3) = noamin
            valk(ik+4) = cpamin
            call tbajli(nomres, nbpar, nopara, vali, valr, &
                        [c16b], valk, 0)
!
        else
            call utmess('F', 'ALGORITH10_56', sk=tych)
        end if
!
100     continue
    end do
!
    call jedetr(knum)
!
999 continue
!
    call jedema()
!
end subroutine
