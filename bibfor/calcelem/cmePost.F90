! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
subroutine cmePost(matrElemz)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/redetr.h"
#include "asterfort/sdmpic.h"
!
character(len=*), intent(in) :: matrElemz
!
! --------------------------------------------------------------------------------------------------
!
! CALC_MATR_ELEM
!
! Post-treatment
!
! --------------------------------------------------------------------------------------------------
!
! In  matrElem         : elementary matrix
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: matrElem, resuElem
    integer :: iexi
    character(len=8) :: mesh, answer
    character(len=24), pointer :: relr(:) => null()
    integer :: iResuElem, nbResuElem
!
! --------------------------------------------------------------------------------------------------
!
    matrElem = matrElemz

! - Get mesh
    call dismoi('NOM_MAILLA', matrElem, 'MATR_ELEM', repk = mesh)

! - Complete for MPI
    call jeexin(matrElem//'.RELR', iexi)
    if (iexi .gt. 0 .and. .not. isParallelMesh(mesh)) then
        call jelira(matrElem//'.RELR', 'LONMAX', nbResuElem)
        call jeveuo(matrElem//'.RELR', 'L', vk24 = relr)
        do iResuElem = 1, nbResuElem
            resuElem = relr(iResuElem)(1:19)
            call jeexin(resuElem//'.RESL', iexi)
            if (iexi .eq. 0) cycle
            call dismoi('MPI_COMPLET', resuElem, 'RESUELEM', repk = answer)
            ASSERT((answer.eq.'OUI').or.(answer.eq.'NON'))
            if (answer .eq. 'NON') then
                call sdmpic('RESUELEM', resuElem)
            endif
        end do
    endif

! - Clean some resuElem 
    call redetr(matrElem)
!
end subroutine
