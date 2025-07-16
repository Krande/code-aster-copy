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
subroutine nmdorc(modelZ, chmateZ, l_etat_init, comporZ, carcriZ, mult_compZ_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/nmdocc.h"
#include "asterfort/nmdocm.h"
#include "asterfort/nmdocr.h"
#include "asterfort/utmess.h"
#include "asterfort/infniv.h"
!
    character(len=*), intent(in) :: modelZ, chmateZ
    aster_logical, intent(in) :: l_etat_init
    character(len=*), intent(in) :: comporZ
    character(len=*), intent(in) :: carcriZ
    character(len=*), optional, intent(in) :: mult_compZ_
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of behaviours (mechanics)
!
! Prepare objects for constitutive laws
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : model
! In  chmate           : material field
! In  l_etat_init      : .true. if initial state is defined
! In  compor           : map for parameters of constitutive laws
! Out carcri           : map for parameters for integration of constitutive law
! In  mult_comp        : map for multi-materials
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=8) :: model, chmate
    character(len=19) :: compor
    character(len=24) :: carcri
!
! --------------------------------------------------------------------------------------------------
!
    model = modelZ
    chmate = chmateZ
    compor = comporZ
    carcri = carcriZ
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE12_2')
    end if

! - Get parameters from COMPORTEMENT keyword and prepare COMPOR map
    call nmdocc(model, chmate, l_etat_init, compor, 'V')

! - Get parameters from COMPORTEMENT keyword and prepare CARCRI map
    call nmdocr(model, carcri, 'V')

! - Get parameters from COMPORTEMENT keyword and prepare MULT_COMP <CARTE> (for crystals)
    if (present(mult_compZ_)) then
        call nmdocm(model, mult_compZ_, 'V')
    end if
!
end subroutine
