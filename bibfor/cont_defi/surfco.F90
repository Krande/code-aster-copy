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
subroutine surfco(sdcont, mesh)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cfdisi.h"
#include "asterfort/infniv.h"
#include "asterfort/surfc1.h"
#include "asterfort/surfc2.h"
#include "asterfort/surfcl.h"
#include "asterfort/surfcp.h"
!
    character(len=8), intent(in) :: sdcont
    character(len=8), intent(in) :: mesh
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Print debug
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  sdcont           : name of contact concept (DEFI_CONTACT)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: unit_msg, lvel_msg
    character(len=24) :: sdcont_defi
    integer(kind=8) :: cont_form
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(unit_msg, lvel_msg)
!
! - Datastructure for contact definition
!
    sdcont_defi = sdcont(1:8)//'.CONTACT'
!
! - Parameters
!
    cont_form = cfdisi(sdcont_defi, 'FORMULATION')
!
! - Debug print
!
    if (lvel_msg .ge. 2) then
        call surfcp(sdcont, unit_msg)
        if (cont_form .eq. 1) then
            call surfcl(sdcont, mesh, unit_msg)
            call surfc1(sdcont, unit_msg)
        else if (cont_form .eq. 2) then
            call surfcl(sdcont, mesh, unit_msg)
            call surfc2(sdcont, mesh)
        end if
    end if
!
end subroutine
