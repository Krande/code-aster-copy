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

subroutine ntcrch(model, nume_dof, vhydr_, hydr_init_)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/alchml.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/vtcreb.h"
!
    character(len=8), intent(in) :: model
    character(len=24), intent(in) :: nume_dof
    character(len=24), optional, intent(in) :: vhydr_
    character(len=24), optional, intent(out) :: hydr_init_
!
! --------------------------------------------------------------------------------------------------
!
! THER_LINEAIRE - Init
!
! Create unknowns
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  nume_dof         : name of numbering object (NUME_DDL)
! In  vhydr            : field for hydration
! Out hydr_init        : field for initial hydration
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: hydric, hydris, ligrmo
    integer(kind=8) :: iret
    character(len=24) :: vtemp
!
! --------------------------------------------------------------------------------------------------
!
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrmo)
!
! - Create hydration
!
    hydric = '&&NTCRCH.HYDR_C'
    hydris = '&&NTCRCH.HYDR_S'
    if (present(vhydr_)) then
        hydr_init_ = '&&NTCRCH.HYDR0'
        call alchml(ligrmo, "TOU_INI_ELGA", "PHYDR_R", "V", hydr_init_, iret, " ")
        ASSERT(iret == 0)
        call copisd('CHAMP_GD', 'V', hydr_init_, vhydr_)
    end if
!
! - Create temperature
!
    vtemp = '&&NXLECTVAR_____'
    call vtcreb(vtemp, 'V', 'R', nume_ddlz=nume_dof)
!
    call detrsd('CHAMP', hydric)
    call detrsd('CHAMP', hydris)
!
end subroutine
