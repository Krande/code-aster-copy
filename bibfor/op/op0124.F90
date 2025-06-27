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

subroutine op0124()
    implicit none
!
!     COMMANDE:  CREA_RESU
!
! ----------------------------------------------------------------------
#include "asterc/getfac.h"
#include "asterfort/crasse.h"
#include "asterfort/crcoch.h"
#include "asterfort/crcore.h"
#include "asterfort/crkucv.h"
#include "asterfort/crperm.h"
#include "asterfort/crprol.h"
#include "asterfort/crtype.h"
#include "asterfort/crvarc.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/ve0124.h"
    integer(kind=8) :: nbfac
    character(len=16) :: typres, valk(2)
!     ------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
!
    call ve0124(typres)
    valk(1) = typres
!
! ----------------------------------------------------------------------
!                   TRAITEMENT DU MOT CLE "PERM_CHAM"
! ----------------------------------------------------------------------
!
    call getfac('PERM_CHAM', nbfac)
    if (nbfac .gt. 0) then
        if (typres .ne. 'EVOL_NOLI') then
            valk(2) = 'PERM_CHAM'
            call utmess('F', 'ALGORITH17_41', nk=2, valk=valk)
        end if
        call crperm()
        goto 999
    end if
!
! ----------------------------------------------------------------------
!               TRAITEMENT DU MOT CLE "PROL_RTZ"
! ----------------------------------------------------------------------
!
    call getfac('PROL_RTZ', nbfac)
    if (nbfac .gt. 0) then
        if (typres .ne. 'EVOL_THER') then
            valk(2) = 'EVOL_THER'
            call utmess('F', 'ALGORITH17_41', nk=2, valk=valk)
        end if
        call crprol()
        goto 999
    end if
!
! ----------------------------------------------------------------------
!               TRAITEMENT DU MOT CLE "AFFE"
! ----------------------------------------------------------------------
!
    call getfac('AFFE', nbfac)
    if (nbfac .gt. 0) then
        call crtype()
        goto 999
    end if
!
! ----------------------------------------------------------------------
!               TRAITEMENT DU MOT CLE "ASSE"
! ----------------------------------------------------------------------
!
    call getfac('ASSE', nbfac)
    if (nbfac .gt. 0) then
        if (typres .ne. 'EVOL_THER') then
            valk(2) = 'ASSE'
            call utmess('F', 'ALGORITH17_41', nk=2, valk=valk)
        end if
        call crasse()
        goto 999
    end if
!
! ----------------------------------------------------------------------
!               TRAITEMENT DU MOT CLE "PREP_VARC"
! ----------------------------------------------------------------------
!
    call getfac('PREP_VARC', nbfac)
    if (nbfac .gt. 0) then
        if (typres .ne. 'EVOL_THER') then
            valk(2) = 'PREP_VARC'
            call utmess('F', 'ALGORITH17_41', nk=2, valk=valk)
        end if
        !
        call crvarc()
        goto 999
    end if
!
! ----------------------------------------------------------------------
!               TRAITEMENT DU MOT CLE "KUCV"
! ----------------------------------------------------------------------
!
    call getfac('KUCV', nbfac)
    if (nbfac .gt. 0) then
        call crkucv()
        goto 999
    end if
!
! ----------------------------------------------------------------------
!               TRAITEMENT DU MOT CLE "CONV_CHAR"
! ----------------------------------------------------------------------
!
    call getfac('CONV_CHAR', nbfac)
    if (nbfac .gt. 0) then
        call crcoch()
        goto 999
    end if
!
!
! ----------------------------------------------------------------------
!               TRAITEMENT DU MOT CLE "CONV_RESU"
! ----------------------------------------------------------------------
!
    call getfac('CONV_RESU', nbfac)
    if (nbfac .gt. 0) then
        call crcore()
        goto 999
    end if
!
999 continue
    call jedema()
end subroutine
