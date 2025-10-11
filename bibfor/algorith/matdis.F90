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

subroutine matdis(matd, verbose)
    implicit none
#include "asterc/getexm.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/utmess.h"
! person_in_charge: jacques.pellet at edf.fr
! ----------------------------------------------------------------------
!
! out K3 matd  : 'OUI' / 'NON'  (MATR_DISTRIBUEE)
!
! ----------------------------------------------------------------------
    character(len=3) :: matd
    aster_logical, intent(in), optional :: verbose
! ----------------------------------------------------------------------
    integer(kind=8) :: ibid, eximc, eximo, iexi
    aster_logical :: verbose_loc
    character(len=8) :: modele, partit

!   -- MATR_DISTRIBUEE ?
    verbose_loc = .true.
    if (present(verbose)) then
        verbose_loc = verbose
    end if
    matd = 'NON'
    eximc = getexm('SOLVEUR', 'MATR_DISTRIBUEE')
    eximo = getexm(' ', 'MODELE')
    if (eximc .eq. 1 .and. eximo .eq. 1) then
        call getvtx('SOLVEUR', 'MATR_DISTRIBUEE', iocc=1, scal=matd, nbret=ibid)
        if (ibid .eq. 0) then
            matd = 'NON'
        end if
        call getvid(' ', 'MODELE', scal=modele, nbret=ibid)
        ASSERT(ibid .eq. 1)
        call dismoi('PARTITION', modele, 'MODELE', repk=partit)
        call exisd('PARTITION', partit, iexi)
        if (iexi .eq. 0 .and. matd .eq. 'OUI') then
            matd = 'NON'
            if (verbose_loc) then
                call utmess('I', 'ASSEMBLA_3')
            end if
        end if
    end if
end subroutine
