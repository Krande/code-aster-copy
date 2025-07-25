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

subroutine op0017()
    implicit none
!
!     COMMANDE:  IMPR_CO
!
!     ------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getltx.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/ulexis.h"
#include "asterfort/ulopen.h"
#include "asterfort/utimsd.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
    integer(kind=8) :: nivo, n3, n1, ifi, n2, nbocc, ncon, ipos, long(1), n4
    integer(kind=8) :: i, iocc
    aster_logical :: lattr, lcont
    character(len=1) :: base
    character(len=3) :: perm
    character(len=8) :: leresu
    character(len=16) :: nomfi
    character(len=72) :: chaine
    character(len=8), pointer :: liste_co(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
!
    nivo = 0
    call getvis(' ', 'NIVEAU', scal=nivo, nbret=n3)
    perm = ''
    if (nivo .eq. -1) call getvtx(' ', 'PERMUTATION', scal=perm)
!
    call getvtx(' ', 'ATTRIBUT', scal=chaine, nbret=n3)
    if (chaine(1:3) .eq. 'OUI') then
        lattr = .true.
    else
        lattr = .false.
    end if
!
    call getvtx(' ', 'CONTENU', scal=chaine, nbret=n3)
    if (chaine(1:3) .eq. 'OUI') then
        lcont = .true.
    else
        lcont = .false.
    end if
!
    call getvtx(' ', 'BASE', scal=chaine, nbret=n1)
    base = chaine(1:1)
!
    ifi = 0
    nomfi = ' '
    call getvis(' ', 'UNITE', scal=ifi, nbret=n2)
    if (.not. ulexis(ifi)) then
        call ulopen(ifi, ' ', nomfi, 'NEW', 'O')
    end if
!
    call getfac('CONCEPT', nbocc)
    do iocc = 1, nbocc
        call getvid('CONCEPT', 'NOM', iocc=iocc, nbval=0, nbret=ncon)
        ncon = -ncon
        if (ncon .gt. 0) then
            AS_ALLOCATE(vk8=liste_co, size=ncon)
            call getvid('CONCEPT', 'NOM', iocc=iocc, nbval=ncon, vect=liste_co, &
                        nbret=n1)
            do i = 1, ncon
                leresu = liste_co(i)
                call utimsd(ifi, nivo, lattr, lcont, leresu, &
                            1, base, perm=perm)
            end do
            AS_DEALLOCATE(vk8=liste_co)
        end if
    end do
!
    call getvtx(' ', 'CHAINE', scal=chaine, nbret=n2)
    if (n2 .gt. 0) then
        call getltx(' ', 'CHAINE', 1, 72, 1, &
                    long, n3)
        call getvis(' ', 'POSITION', scal=ipos, nbret=n4)
        call utimsd(ifi, nivo, lattr, lcont, chaine(1:long(1)), &
                    ipos, base, perm=perm)
    end if
!
    call getvtx(' ', 'TOUT', scal=chaine, nbret=n2)
    if (n2 .gt. 0) then
        call utimsd(ifi, nivo, lattr, lcont, ' ', &
                    0, base, perm=perm)
    end if
    flush (ifi)
    if (ifi .ne. 6) then
        call ulopen(-ifi, ' ', ' ', ' ', ' ')
    end if
    call jedema()
end subroutine
