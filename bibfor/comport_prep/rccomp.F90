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

subroutine rccomp(chmat, mesh)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/comp_comp_read.h"
#include "asterfort/comp_comp_save.h"
#include "asterfort/comp_init.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
!
    character(len=8), intent(in) :: chmat
    character(len=8), intent(in) :: mesh
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (AFFE_MATERIAU)
!
! Prepare COMPOR <CARTE>
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh        : name of mesh
! In  chmat       : name material field
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nocc
    character(len=16), parameter :: factorKeyword = 'AFFE_COMPOR'
    character(len=19) :: compor
    character(len=16), pointer :: infoValk(:) => null()
    integer(kind=8), pointer :: infoVali(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    compor = chmat//'.COMPOR'
    call getfac(factorKeyword, nocc)
    if (nocc .ne. 0) then

! ----- Create comportment informations objects
        AS_ALLOCATE(vk16=infoValk, size=16*nocc)
        AS_ALLOCATE(vi=infoVali, size=4*nocc)

! ----- Create COMPOR <CARTE>
        call comp_init(mesh, compor, 'G')

! ----- Read informations from command file
        call comp_comp_read(infoValk, infoVali)

! ----- Save informations in COMPOR <CARTE>
        call comp_comp_save(mesh, compor, infoValk, infoVali)

! ----- Clean
        AS_DEALLOCATE(vk16=infoValk)
        AS_DEALLOCATE(vi=infoVali)
    end if
!
end subroutine
