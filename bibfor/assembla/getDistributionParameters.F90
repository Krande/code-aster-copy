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
subroutine getDistributionParameters(nbElem, listElem, &
                                     ldist, ldgrel, &
                                     rang, nbproc, &
                                     numsd)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/jeveuo.h"
#include "asterfort/parti0.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: nbElem
    character(len=*), intent(in) :: listElem(nbElem)
    aster_logical, intent(out) :: ldist, ldgrel
    integer(kind=8), intent(out) :: rang, nbproc
    integer(kind=8), pointer :: numsd(:)
!
! --------------------------------------------------------------------------------------------------
!
    !   Get parameters for distribution of elementary  vectors/matrices
!
! --------------------------------------------------------------------------------------------------
!
    mpi_int :: mrank, msize
    character(len=8) :: partit
    character(len=24), pointer :: prtk(:) => null()
    integer(kind=8), pointer :: prti(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    numsd => null()
    ldist = ASTER_FALSE
    ldgrel = ASTER_FALSE
    rang = 0
    nbproc = 1

! - Get partition parameter for elementary terms
    call parti0(nbElem, listElem, partit)

    if (partit .ne. ' ') then
        ldist = ASTER_TRUE
        call asmpi_info(rank=mrank, size=msize)
        rang = to_aster_int(mrank)
        nbproc = to_aster_int(msize)
        call jeveuo(partit//'.PRTK', 'L', vk24=prtk)
        ldgrel = prtk(1) .eq. 'SOUS_DOMAINE' .or. prtk(1) .eq. 'GROUP_ELEM'
        if (.not. ldgrel) then
            call jeveuo(partit//'.PRTI', 'L', vi=prti)
            if (prti(1) .gt. nbproc) then
                call utmess('F', 'CALCUL_35', ni=2, vali=[prti(1), nbproc])
            end if
            call jeveuo(partit//'.NUPR', 'L', vi=numsd)
        end if
    end if

end subroutine
