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
subroutine check_model(mesh, cont_form)
!
    implicit none
!
#include "asterfort/exipat.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(in) :: cont_form
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Check if  PATCH in mesh
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  cont_form        : formulation of contact
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret
!
! --------------------------------------------------------------------------------------------------
!

! - Check if exist PATCH in mesh (LAC method)
    if (cont_form .eq. 5) then
        call exipat(mesh, iret)
        if (iret .eq. 0) then
            call utmess('F', 'CONTACT4_2', sk=mesh)
        end if
    end if
!
end subroutine
