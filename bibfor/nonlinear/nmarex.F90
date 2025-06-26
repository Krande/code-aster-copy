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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine nmarex(keywfact, sdarch, lDyna_)
!
    implicit none
!
#include "asterfort/getvtx.h"
#include "asterfort/wkvect.h"
#include "asterfort/utmess.h"
!
    character(len=19), intent(in) :: sdarch
    character(len=16), intent(in) :: keywfact
    aster_logical, optional, intent(in) :: lDyna_
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE *_NON_LINE (ARCHIVAGE)
!
! CONSTRUCTION CHAMPS EXCLUS DE L'ARCHIVAGE
!
! --------------------------------------------------------------------------------------------------
!
! IN  MOTFAC : MOT-FACTEUR POUR LIRE <CHAM_EXCL>
! IN  SDARCH : NOM DE LA SD ARCHIVAGE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbFieldExcl, iFieldExcl
    aster_logical :: lDyna, lAcce, lSigm
    character(len=16), pointer :: fieldExcl(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    lDyna = ASTER_FALSE
    if (present(lDyna_)) then
        lDyna = lDyna_
    end if
!
    call getvtx(keywfact, 'CHAM_EXCLU', iocc=1, nbval=0, nbret=nbFieldExcl)
    nbFieldExcl = -nbFieldExcl
    if (nbFieldExcl .ne. 0) then
        call wkvect(sdarch(1:19)//'.AEXC', 'V V K16', nbFieldExcl, vk16=fieldExcl)
        call getvtx(keywfact, 'CHAM_EXCLU', iocc=1, nbval=nbFieldExcl, vect=fieldExcl)
    end if
!
    if (lDyna) then
        lSigm = ASTER_TRUE
        lAcce = ASTER_TRUE
        do iFieldExcl = 1, nbFieldExcl
            if (fieldExcl(iFieldExcl) .eq. 'SIEF_ELGA') then
                lSigm = ASTER_FALSE
            end if
            if (fieldExcl(iFieldExcl) .eq. 'ACCE') then
                lAcce = ASTER_FALSE
            end if
        end do
        if (lSigm .and. .not. lAcce) then
            call utmess('A', 'MECANONLINE5_26')
        end if
    end if

!
end subroutine
