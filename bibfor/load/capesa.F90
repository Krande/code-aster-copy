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

subroutine capesa(load, mesh, valeType, nbOcc)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterfort/assert.h"
#include "asterfort/getelem.h"
#include "asterfort/char_crea_cart.h"
#include "asterfort/getvr8.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/nocart.h"
!
    character(len=8), intent(in) :: load, mesh
    character(len=4), intent(in) :: valeType
    integer(kind=8), intent(in) :: nbOcc
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Keyword = 'PESANTEUR'
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  mesh             : mesh
! In  valeType         : affected value type (real, complex or function)
! In  nbOcc            : number of factor keywords
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordfact = 'PESANTEUR'
    character(len=24), parameter :: listCell = '&&CAPESA.LIST_ELEM'
    real(kind=8) :: pesa(4), norme, pes(3)
    integer(kind=8) :: iocc, npesa
    integer(kind=8) :: jvCell, nbCell
    real(kind=8), pointer :: valv(:) => null()
    character(len=19) :: map(LOAD_MAP_NBMAX)
    integer(kind=8) :: nbMap, nbCmp(LOAD_MAP_NBMAX)
!
! --------------------------------------------------------------------------------------------------
!

!
! - Initializations
!
    ASSERT(valeType .eq. 'REEL')
!
! - Creation and initialization to zero of <CARTE>
!
    call char_crea_cart('MECANIQUE', keywordfact, load, mesh, valeType, &
                        nbMap, map, nbCmp)
    ASSERT(nbMap .eq. 1)
    call jeveuo(map(1)//'.VALV', 'E', vr=valv)
!
! - Loop on factor keyword
!
    do iocc = 1, nbOcc

! ----- Get values of load
        call getvr8('PESANTEUR', 'GRAVITE', iocc=iocc, scal=pesa(1), nbret=npesa)
        call getvr8('PESANTEUR', 'DIRECTION', iocc=iocc, nbval=3, vect=pes, nbret=npesa)
        norme = sqrt(pes(1)*pes(1)+pes(2)*pes(2)+pes(3)*pes(3))
        if (norme .gt. r8miem()) then
            pesa(2) = pes(1)/norme
            pesa(3) = pes(2)/norme
            pesa(4) = pes(3)/norme
        else
            call utmess('F', 'CHARGES2_53')
        end if
        valv(1:4) = pesa(1:4)

! ----- Read mesh and affect
        call getelem(mesh, keywordfact, iocc, ' ', listCell, nbCell)
        if (nbCell .eq. 0) then
            call nocart(map(1), 1, nbCmp(1))
        else
            call jeveuo(listCell, 'L', jvCell)
            call nocart(map(1), 3, nbCmp(1), mode='NUM', nma=nbCell, &
                        limanu=zi(jvCell))
        end if

    end do
end subroutine
