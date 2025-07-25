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
subroutine rvvsup()
    implicit none
!     VERIFICATION SUPPLEMENTAIRE OP0051
!     ------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/gettco.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: n1, n2, n3, n4, iocc, nbpost
    character(len=8) :: resu, nomres
    character(len=24) :: valk(4)
    character(len=16) :: nomcmd, concep, typres
!
!=======================================================================
!
    call jemarq()
    call getres(resu, concep, nomcmd)
!
!     --- VERIFICATION SUR "OPERATION" ---
!
    call getfac('ACTION', nbpost)
!
    do iocc = 1, nbpost, 1
!
!     /* QUANTITE (IE : SOMME) */
        call getvtx('ACTION', 'RESULTANTE', iocc=iocc, nbval=0, nbret=n1)
        n1 = -n1
        if (n1 .gt. 0) then
            call getvtx('ACTION', 'RESULTANTE', iocc=iocc, nbval=0, nbret=n1)
            call getvtx('ACTION', 'MOMENT    ', iocc=iocc, nbval=0, nbret=n2)
            n1 = -n1
            n2 = -n2
            if (n2 .ne. 0) then
                if (((n1 .ne. 2) .and. (n1 .ne. 3)) .or. (n1 .ne. n2)) then
                    call utmess('F', 'POSTRELE_42', si=iocc)
                end if
                call getvr8('ACTION', 'POINT', iocc=iocc, nbval=0, nbret=n1)
                n1 = -n1
                if ((n1 .ne. 2) .and. (n1 .ne. 3)) then
                    call utmess('F', 'POSTRELE_43', si=iocc)
                end if
            end if
        end if
!
!     /* COHERENCE ACCES DANS RESULTAT */
        call getvid('ACTION', 'RESULTAT', iocc=iocc, nbval=0, nbret=n1)
        n1 = -n1
        if (n1 .gt. 0) then
            call getvid('ACTION', 'RESULTAT', iocc=iocc, scal=nomres, nbret=n1)
            call gettco(nomres, typres)
            call getvid('ACTION', 'LIST_FREQ', iocc=iocc, nbval=0, nbret=n1)
            call getvr8('ACTION', 'FREQ', iocc=iocc, nbval=0, nbret=n2)
            n1 = max(-n1, -n2)
            call getvid('ACTION', 'LIST_INST', iocc=iocc, nbval=0, nbret=n2)
            call getvr8('ACTION', 'INST', iocc=iocc, nbval=0, nbret=n3)
            n2 = max(-n3, -n2)
            call getvid('ACTION', 'LIST_MODE', iocc=iocc, nbval=0, nbret=n3)
            call getvis('ACTION', 'NUME_MODE', iocc=iocc, nbval=0, nbret=n4)
            n3 = max(-n3, -n4)
            n4 = max(n1, n2, n3)
            if (n4 .gt. 0) then
                if (((n1 .ne. 0) .or. (n3 .ne. 0)) .and. &
                    ( &
                    (typres(1:4) .eq. 'EVOL') .or. (typres(6:10) .eq. 'TRANS') .or. &
                    (typres(11:15) .eq. 'TRANS') &
                    )) then
                    valk(1) = nomres
                    valk(2) = typres
                    valk(3) = 'FREQ'
                    valk(4) = 'MODE'
                    call utmess('F', 'POSTRELE_44', nk=4, valk=valk, si=iocc)
                end if
                if ((n2 .ne. 0) .and. &
                    ( &
                    (typres(1:4) .eq. 'MODE') .or. (typres(1:4) .eq. 'BASE') .or. &
                    (typres(6:10) .eq. 'HARMO') .or. (typres(11:15) .eq. 'HARMO') &
                    )) then
                    valk(1) = nomres
                    valk(2) = typres
                    valk(3) = 'INSTANT'
                    call utmess('F', 'POSTRELE_45', nk=3, valk=valk, si=iocc)
                end if
            end if
        end if
!
    end do
!
    call jedema()
end subroutine
