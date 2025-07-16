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

subroutine ibtcpu(ier)
    implicit none
#include "asterc/gtoptr.h"
#include "asterc/jdcget.h"
#include "asterc/rdtmax.h"
#include "asterfort/assert.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ier
!     OPTION DE MODIFICATION DE LA LIMITE DE TEMPS CPU POUR CONSERVER
!     UNE MARGE SUFFISANTE AFIN DE TERMINER PROPREMENT UNE EXECUTION
!     ------------------------------------------------------------------
!            0 TOUT C'EST BIEN PASSE
!            1 ERREUR DANS LA LECTURE DE LA COMMANDE
!     ------------------------------------------------------------------
!
    integer(kind=8) :: l1, l2, lcpu, iret, vali(3), itest
    real(kind=8) :: pccpu, tpmax, dix, ntmax
    parameter(dix=10.d0)
!
    ier = 0
    tpmax = 0.d0
    ntmax = 0.d0
    l1 = 0
    l2 = 0
!     RECUPERATION DU TEMPS LIMITE DE L'EXECUTION
    call gtoptr('tpmax', tpmax, iret)
    ASSERT(iret .eq. 0)
!
    itest = jdcget('TestMode')
!
    call getvis('RESERVE_CPU', 'VALE', iocc=1, scal=lcpu, nbret=l1)
    call getvr8('RESERVE_CPU', 'POURCENTAGE', iocc=1, scal=pccpu, nbret=l2)
!
    if (l1 .eq. 0 .and. l2 .eq. 0) then
!       set default values here
        if (itest .eq. 0) then
            ntmax = tpmax*0.9d0
        else
            ntmax = tpmax-dix
        end if
    elseif (l1 .ne. 0) then
        if (lcpu .gt. tpmax) then
            call utmess('F', 'SUPERVIS_31')
            ier = 1
        end if
        ntmax = tpmax-lcpu
    else
!       l2 .ne. 0
        ntmax = tpmax*(1.0-pccpu)
    end if

!   reduce execution time limit
    call rdtmax(tpmax-ntmax)
    vali(1) = int(tpmax)
    vali(2) = int(ntmax)
    vali(3) = int(tpmax-ntmax)
    call utmess('I', 'SUPERVIS_64', ni=3, vali=vali)
!
end subroutine
