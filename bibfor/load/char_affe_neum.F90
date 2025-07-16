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
subroutine char_affe_neum(model, mesh, ndim, &
                          keywordfact, iocc, &
                          nbMap, map, nbCmp)
!
    implicit none
!
#include "asterfort/getelem.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/vetyma.h"
!
    character(len=8), intent(in) :: model, mesh
    integer(kind=8), intent(in) :: ndim
    character(len=16), intent(in) :: keywordfact
    integer(kind=8), intent(in) :: iocc, nbMap
    character(len=19), intent(in) :: map(nbMap)
    integer(kind=8), intent(in) :: nbCmp(nbMap)
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Apply Neumann loads in <CARTE> with elements type check
!
! --------------------------------------------------------------------------------------------------
!
! In  model        : model
! In  mesh         : mesh
! In  ndim         : space dimension
! In  keywordfact  : factor keyword to read elements
! In  iocc         : factor keyword index in AFFE_CHAR_MECA
! In  nb_carte     : number of <CARTE> for this Neumann load
! In  carte        : <CARTE> for this Neumann load
! In  nb_cmp       : number of components in the <CARTE>
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24), parameter:: listCell = '&&LIST_ELEM'
    integer(kind=8), pointer :: cellNume(:) => null()
    integer(kind=8) :: nbCell
    integer(kind=8) :: iMap
!
! --------------------------------------------------------------------------------------------------
!

!
! - Get list of cells
!
    call getelem(mesh, keywordfact, iocc, 'A', listCell, nbCell, model=model)

    if (nbCell .ne. 0) then

! ----- Check elements
        call jeveuo(listCell, 'L', vi=cellNume)
        do iMap = 1, nbMap
            if (nbCmp(iMap) .ne. 0) then
                call vetyma(mesh, ndim, keywordfact, listCell, nbCell)
            end if
        end do

! ----- Apply Neumann loads in <CARTE>
        do iMap = 1, nbMap
            if (nbCmp(iMap) .ne. 0) then
                call nocart(map(iMap), 3, nbCmp(iMap), mode='NUM', nma=nbCell, &
                            limanu=cellNume)
            end if
        end do
!
    end if
!
    call jedetr(listCell)
!
end subroutine
