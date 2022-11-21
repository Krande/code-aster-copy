! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine updateMateco(mater, mateco, l_thm, l_ther)
!
implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/jedetc.h"
!
character(len=8), intent(in) :: mateco, mater
aster_logical, intent(in) :: l_thm, l_ther
!
! --------------------------------------------------------------------------------------------------
!
! Material
!
! Update jeveux pointers if codedMaterial (called from c++)
!
! --------------------------------------------------------------------------------------------------
!
! In  mateco           : name of coded material (MATER_CODE)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: k19bid
!
! --------------------------------------------------------------------------------------------------
!
!
    call jedetc(' ', mateco, 1)
!
    call rcmfmc(mater, k19bid, l_thm, l_ther, mateco, 'G')
!
end subroutine
