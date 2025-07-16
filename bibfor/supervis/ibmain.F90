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

subroutine ibmain()
    use logging_module, only: initialize
    use parameters_module, only: ST_OK
    use superv_module, only: superv_before
    implicit none

#include "asterc/gtopti.h"
#include "asterc/gtoptr.h"
#include "asterc/inisig.h"
#include "asterc/ismaem.h"
#include "asterc/loisem.h"
#include "asterfort/entete.h"
#include "asterfort/ibimpr.h"
#include "asterfort/jedebu.h"
#include "asterfort/jeinif.h"
#include "asterfort/lxinit.h"
#include "asterfort/ststat.h"
#include "asterfort/utmess.h"

    character(len=8) :: nomf
    integer(kind=8) :: unmega, idebug, iret, lois
    integer(kind=8) :: mxdyn
    real(kind=8) :: valr(2), moctet, memory

!   Initialization of loggers
    call initialize()

!   Initialization of the internal parser
    call lxinit()

!   Initialization of signal interruption
    call inisig()

!   Initialization of the global status
    call ststat(ST_OK)

!   Initialization of logical units
    call ibimpr()

!   Initialization of jeveux
    idebug = 0
    call gtopti('dbgjeveux', idebug, iret)
    memory = 0.d0
    call gtoptr('memory', memory, iret)
    unmega = 1024*1024
    lois = loisem()
    moctet = memory*unmega
    valr = 0.d0
    if (moctet .gt. ismaem()) then
        valr(1) = moctet
        valr(2) = ismaem()
        call utmess('F', 'JEVEUX_1', nr=2, valr=valr)
    end if
    mxdyn = int(moctet)
    call jedebu(4, mxdyn/lois, idebug)

!   Allocate a temporary volatile database
    nomf = 'VOLATILE'
    call jeinif('DEBUT', 'DETRUIT', nomf, 'V', 250, 100, 1)
    call superv_before()

!   Print header, without the memory informations
!   (it will be done by adjust_memlimit after ibcata)
    call entete()

    if (idebug .eq. 1) then
        call utmess('I', 'SUPERVIS_12')
    end if

end subroutine
