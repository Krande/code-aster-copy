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

subroutine mdlibe(nomres, nbnli, sd_index)
    implicit none
!
! person_in_charge: hassan.berro at edf.fr
!
!     Memory clearup of resulting vectors in a DYNA_VIBRA//TRAN/GENE calculation
!     with a non-constant integration step
!     ------------------------------------------------------------------
! in  : nomres : result name (usually &&AD****)
! in  : nbnli  : number of localised non linearities
! ----------------------------------------------------------------------
!
#include "asterfort/codent.h"
#include "asterfort/jelibe.h"
#include "asterfort/jeexin.h"
!   Input arguments
    character(len=8), intent(in) :: nomres
    integer(kind=8), intent(in) :: nbnli
    integer(kind=8), optional, intent(in) :: sd_index
!   Local variables
    integer(kind=8)           :: iret, sd_ind
    character(len=4)  :: bl3pt
    character(len=7)  :: intk7
    character(len=16) :: nomres16
!-----------------------------------------------------------------------
    bl3pt = '   .'
    sd_ind = 0

    if (present(sd_index)) sd_ind = sd_index

!   Hard-encode the sd index into the result name
    nomres16 = nomres//'        '
    if (sd_ind .gt. 0) then
        call codent(sd_ind, 'D0', intk7)
        nomres16 = nomres//'.'//intk7
    end if

    call jelibe(nomres16//bl3pt//'DEPL')
    call jelibe(nomres16//bl3pt//'VITE')
    call jelibe(nomres16//bl3pt//'ACCE')
    call jelibe(nomres16//bl3pt//'ORDR')
    call jelibe(nomres16//bl3pt//'DISC')
    call jelibe(nomres16//bl3pt//'PTEM')

    if (nbnli .gt. 0) then
        call jeexin(nomres16//'.NL.VINT', iret)
        if (iret .gt. 0) call jelibe(nomres16//'.NL.VINT')
    end if
!
end subroutine
