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
subroutine nxdocc(model, compor, base_)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/comp_init.h"
#include "asterfort/comp_ther_read.h"
#include "asterfort/comp_ther_save.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedetr.h"
!
    character(len=8), intent(in) :: model
    character(len=19), intent(in) :: compor
    character(len=1), optional, intent(in) :: base_
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (thermics)
!
! Get parameters from COMPORTEMENT keyword and prepare COMPOR <CARTE>
!
! --------------------------------------------------------------------------------------------------
!
! In  model       : name of model
! In  compor      : name of <CARTE> COMPOR
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbCmp
    character(len=8) :: mesh
    character(len=1) :: base
    character(len=19), parameter :: list_vale = '&&LIST_VALE'
!
! --------------------------------------------------------------------------------------------------
!
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
!
    base = "V"
    if (present(base_)) then
        base = base_
    end if
!
! - Read informations from command file
    call comp_ther_read(list_vale)

! - Create COMPOR <CARTE>
    call comp_init(mesh, compor, base, nbCmp)

! - Save informations in COMPOR <CARTE>
    call comp_ther_save(mesh, compor, nbCmp, list_vale)

! - Clean
    call jedetr(list_vale)
!
end subroutine
