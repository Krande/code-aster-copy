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
subroutine nmsssv(modelz, matez, caraElemz, listLoad, vesstf)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/vemare.h"
#include "asterfort/ss2mme.h"
!
    character(len=*), intent(in) :: modelz, matez, caraElemz
    character(len=19), intent(in) :: vesstf, listLoad
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (CALCUL - SOUS-STRUCTURATION)
!
! CALCUL DU VECTEUR CHARGEMENT SUR MACRO-ELEMENTS
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1), parameter :: base = 'V'
    character(len=8) :: model, mate
    character(len=24) :: funcMultSuper, caraElem
    integer(kind=8) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    mate = matez
    caraElem = caraElemz
    model = modelz
    funcMultSuper = listLoad(1:19)//'.FCSS'

! - CALCUL
    call jeexin(funcMultSuper, iret)
    ASSERT(iret .ne. 0)
    call vemare(base, vesstf, model)
    call jedetr(vesstf//'.RELC')
    call ss2mme(model, vesstf, base)
!
    call jedema()
!
end subroutine
