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
subroutine nmdocr(model, carcri, base)
!
    use BehaviourPrepare_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/carc_chck.h"
#include "asterfort/carc_delete.h"
#include "asterfort/carc_info.h"
#include "asterfort/carc_init.h"
#include "asterfort/carc_read.h"
#include "asterfort/carc_save.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in)  :: model
    character(len=24), intent(in) :: carcri
    character(len=1), intent(in) :: base
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of behaviours (mechanics)
!
! Get parameters from COMPORTEMENT keyword and prepare CARCRI <CARTE>
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : model
! In  carcri           : map for parameters for integration of constitutive law
! In  base             : Jeveux base
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=8) :: mesh
    type(BehaviourPrep_MapCarcri) :: prepMapCarcri
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE12_5')
    end if

! - Initialisations
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)

! - Create objects
    call carc_info(prepMapCarcri)

! - Create CARCRI <CARTE>
    call carc_init(mesh, carcri, base)

! - Read informations from command file
    call carc_read(prepMapCarcri, model)

! - Check informations in CARCRI <CARTE>
    call carc_chck(prepMapCarcri)

! - Save and check informations in CARCRI <CARTE>
    call carc_save(mesh, carcri, prepMapCarcri)

! - Cleaning
    call carc_delete(prepMapCarcri)

!
end subroutine
