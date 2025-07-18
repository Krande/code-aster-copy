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
subroutine dltp0(t0, nume)
    implicit none
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rs_getlast.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsorac.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    real(kind=8) :: t0
    integer(kind=8) :: nume
! OUT : T0   : INSTANT INITIAL
! OUT : NUME : NUMERO D'ORDRE DE REPRISE
!     ------------------------------------------------------------------
    integer(kind=8) :: vali
    real(kind=8) :: valr
    character(len=8) :: k8b, nomres, dyna, li, crit, ctype
    character(len=16) :: typres, nomcmd
    character(len=24) :: valk
    complex(kind=8) :: c16b
!     -----------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibid, jadr, jordr, n1, nbordr(1), tnume(1)
    integer(kind=8) :: nbtrou, nc, ndy, nni, np, nt
    real(kind=8) :: prec, r8b, temps
    real(kind=8), pointer :: bint(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    call getres(nomres, typres, nomcmd)
!
!     --- EST-ON EN REPRISE ? ---
!
!     INITIALISATION DE T0 PAR DEFAUT
!
    call getvid('ETAT_INIT', 'RESULTAT', iocc=1, scal=dyna, nbret=ndy)
    if (ndy .ne. 0) then
        call getvis('ETAT_INIT', 'NUME_ORDRE', iocc=1, scal=nume, nbret=nni)
        if (nni .eq. 0) then
            call getvr8('ETAT_INIT', 'INST_INIT', iocc=1, scal=temps, nbret=nt)
            if (nt .eq. 0) then
                call rs_getlast(dyna, nume)
            else
                call getvr8('ETAT_INIT', 'PRECISION', iocc=1, scal=prec, nbret=np)
                call getvtx('ETAT_INIT', 'CRITERE', iocc=1, scal=crit, nbret=nc)
                call rsorac(dyna, 'INST', ibid, temps, k8b, &
                            c16b, prec, crit, tnume, 1, &
                            nbtrou)
                nume = tnume(1)
                if (nbtrou .lt. 0) then
                    valk = dyna
                    valr = temps
                    vali = -nbtrou
                    call utmess('F', 'ALGORITH12_83', sk=valk, si=vali, sr=valr)
                else if (nbtrou .eq. 0) then
                    valk = dyna
                    valr = temps
                    call utmess('F', 'ALGORITH12_84', sk=valk, sr=valr)
                end if
            end if
        else
!           --- VERIFICATION QUE NUME EXISTE ---
            call rsorac(dyna, 'LONUTI', 0, r8b, k8b, &
                        c16b, r8b, k8b, nbordr, 1, &
                        ibid)
            call wkvect('&&COMDLT.NUME_ORDRE', 'V V I', nbordr(1), jordr)
            call rsorac(dyna, 'TOUT_ORDRE', 0, r8b, k8b, &
                        c16b, r8b, k8b, zi(jordr), nbordr(1), &
                        ibid)
            do i = 1, nbordr(1)
                if (zi(jordr+i-1) .eq. nume) goto 12
            end do
            call utmess('F', 'ALGORITH3_36', sk=dyna)
12          continue
        end if
!
!        --- RECUPERATION DE L'INSTANT ---
        call rsadpa(dyna, 'L', 1, 'INST', nume, &
                    1, sjv=jadr, styp=ctype)
        t0 = zr(jadr)
    else
!
!     --- DEFINITION DES INSTANTS DE CALCUL A PARTIR DE "LIST_INST" ---
!
        call getvid('INCREMENT', 'LIST_INST', iocc=1, scal=li, nbret=n1)
        if (n1 .ne. 0) then
            call jeveuo(li//'           .BINT', 'L', vr=bint)
            t0 = bint(1)
        else
!
!
!     --- DEFINITION DE L'INSTANT INITIAL AVEC "INST_INIT" ---
!
            t0 = 0.d0
            call getvr8('INCREMENT', 'INST_INIT', iocc=1, scal=t0, nbret=np)
        end if
    end if
!
    call jedema()
end subroutine
