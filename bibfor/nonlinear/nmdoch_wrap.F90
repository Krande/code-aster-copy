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
subroutine nmdoch_wrap(listLoadZ, jvBaseZ)
!
    use listLoad_type
    use listLoad_module
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/nmdoch.h"
!
    character(len=*), intent(in) :: listLoadZ
    character(len=*), intent(in) :: jvBaseZ
!
! --------------------------------------------------------------------------------------------------
!
! Mechanics - Read parameters
!
! Get loads information and create datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  listLoad          : name of datastructure for list of loads
! In  jvBase            : JEVEUX base where to create objects
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1) :: jvBase
    character(len=24) :: listLoad
    type(ListLoad_Prep) :: listLoadPrep
!
! --------------------------------------------------------------------------------------------------
!
    listLoad = listLoadZ
    jvBase = jvBaseZ

! - Standard loads
    call nmdoch(listLoadPrep, listLoad, jvBase)

! - Debug
    !call listLoadDebug(listLoad)
!
end subroutine
