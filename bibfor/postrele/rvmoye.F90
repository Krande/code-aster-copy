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

subroutine rvmoye(nomres, iocc)
!
!     COMMANDE : POST_RELEVE, OPERATION='MOYENNE_ARITH'
!
! ----------------------------------------------------------------------
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
#include "asterfort/prmono.h"
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
    integer(kind=8) :: nbpar
    character(len=16) :: nopara(200)
!
    integer(kind=8) :: ibid, n1, np, nc, iret, icmp, nbcmp
    integer(kind=8) :: jordr, i100, nbordr, iord, vali(20), nbc
    integer(kind=8) :: ii, ik, ir, jaces, nbacc
    integer(kind=8) :: iadr, iac
    real(kind=8) :: prec, som(64), valr(200)
    complex(kind=8) :: c16b
    character(len=3) :: typpar
    character(len=8) :: crit, resu, nocmp(64), tych, ctype
    character(len=16) :: nomcha, intitu
    character(len=19) :: knum, champ
    character(len=24) :: nomjv
    character(len=80) :: valk(200)
    aster_logical :: exist
!
! ---------------------------------------------------------------------
!
    call jemarq()
    c16b = (0.d0, 0.d0)
    knum = '&&RVMOYE.NUME_ORDRE'
    nbc = 0
!
    call getvtx('ACTION', 'INTITULE', iocc=iocc, scal=intitu, nbret=n1)
    nbpar = 1
    nopara(nbpar) = 'INTITULE'
    valk(1) = intitu
!
! ----- TRAITEMENT DU CHAMP_GD  -----
!
    call getvid('ACTION', 'CHAM_GD', iocc=iocc, scal=champ, nbret=n1)
    if (n1 .ne. 0) then
        nbpar = nbpar+1
        nopara(nbpar) = 'CHAM_GD'
        valk(2) = champ
        call dismoi('TYPE_CHAMP', champ, 'CHAMP', repk=tych)
!
        if (tych(1:4) .eq. 'NOEU') then
            call prmono(champ, iocc, som, nbcmp, nocmp)
            nbpar = nbpar+1
            nopara(nbpar) = 'CMP'
            nbpar = nbpar+1
            nopara(nbpar) = 'MOYENNE'
            do icmp = 1, nbcmp
                valk(3) = nocmp(icmp)
                valr(1) = som(icmp)
                call tbajli(nomres, nbpar, nopara, vali, valr, &
                            [c16b], valk, 0)
            end do
!
        else if (tych(1:2) .eq. 'EL') then
            call utmess('F', 'ALGORITH17_5')
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
    nbpar = nbpar+1
    nopara(nbpar) = 'RESU'
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
    nbpar = nbpar+1
    nopara(nbpar) = 'NOM_CHAM'
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
        nomjv = '&&RVMOYE.NOMS_ACCES'
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
!
        call rsexch(' ', resu, nomcha, iord, champ, &
                    iret)
        if (iret .ne. 0) goto 101
        call dismoi('TYPE_CHAMP', champ, 'CHAMP', repk=tych)
!
        if (tych(1:4) .eq. 'NOEU') then
!
            call prmono(champ, iocc, som, nbcmp, nocmp)
!
            do icmp = 1, nbcmp
                nbpar = nbpar+1
                nopara(nbpar) = 'CMP'
                ik = ik+1
                valk(ik) = nocmp(icmp)
                nbpar = nbpar+1
                nopara(nbpar) = 'MOYENNE'
                ir = ir+1
                valr(ir) = som(icmp)
                call tbajli(nomres, nbpar, nopara, vali, valr, &
                            [c16b], valk, 0)
            end do
!
        else if (tych(1:2) .eq. 'EL') then
            call utmess('F', 'ALGORITH17_5')
        else
            call utmess('F', 'ALGORITH10_56', sk=tych)
        end if
!
101     continue
    end do
!
    call jedetr(knum)
!
999 continue
!
    call jedema()
!
end subroutine
