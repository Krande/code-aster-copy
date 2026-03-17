! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine listun(mesh, zoneKeyword, nbUnilZone, &
                  noponoJv, nolinoJv, nbNodeUnil)
!
    implicit none
!
#include "asterfort/exnode.h"
#include "asterfort/nbnode.h"
!
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(in) :: nbUnilZone
    character(len=16), intent(in) :: zoneKeyword
    character(len=24), intent(in) :: noponoJv, nolinoJv
    integer(kind=8), intent(out) :: nbNodeUnil
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Get nodes from definition of LIAISON_UNILATER
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  nbUnilZone       : number of zones with LIAISON_UNILATER
! Out nbNodeUnil       : total number of nodes in LIAISON_UNILATER
!
! --------------------------------------------------------------------------------------------------
!
    call nbnode(mesh, zoneKeyword, nbUnilZone, noponoJv, nbNodeUnil)
    call exnode(mesh, zoneKeyword, nbUnilZone, nbNodeUnil, nolinoJv)
!
end subroutine
