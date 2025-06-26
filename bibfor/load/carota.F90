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
subroutine carota(load, mesh, valeType)
!
    implicit none
!
#include "asterf_types.h"
#include "LoadTypes_type.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/r8miem.h"
#include "asterfort/assert.h"
#include "asterfort/char_crea_cart.h"
#include "asterfort/char_read_val.h"
#include "asterfort/char_read_vect.h"
#include "asterfort/getelem.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/normev.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: load, mesh
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Keyword = 'ROTATION'
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh      : name of mesh
! In  load      : name of load
! In  valeType  : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordfact = 'ROTATION'
    character(len=24), parameter :: listCell = '&&CAROTA.LIST_ELEM'
    complex(kind=8) :: c16dummy
    character(len=8) :: k8dummy
    character(len=16) :: k16dummy
    real(kind=8) :: rota_speed, rota_axis(3), rota_cent(3), norme
    integer(kind=8) :: iocc, nrota, val_nb
    integer(kind=8) :: jvCell, nbCell
    real(kind=8), pointer :: valv(:) => null()
    character(len=19) :: map(LOAD_MAP_NBMAX)
    integer(kind=8) :: nbMap, nbCmp(LOAD_MAP_NBMAX)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call getfac(keywordfact, nrota)
    if (nrota .eq. 0) goto 99
    ASSERT(nrota .eq. 1)
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
! - Loop on keywords
!
    do iocc = 1, nrota

! ----- Get speed
        call char_read_val(keywordfact, iocc, 'VITESSE', valeType, val_nb, &
                           rota_speed, k8dummy, c16dummy, k16dummy)
        ASSERT(val_nb .eq. 1)

! ----- Get axis
        call char_read_vect(keywordfact, iocc, 'AXE', rota_axis)
        call normev(rota_axis, norme)
        if (norme .le. r8miem()) then
            call utmess('F', 'CHARGES2_53')
        end if

! ----- Get center
        call char_read_vect(keywordfact, iocc, 'CENTRE', rota_cent)

! ----- Affectation of values in <CARTE>
        valv(1) = rota_speed
        valv(2) = rota_axis(1)
        valv(3) = rota_axis(2)
        valv(4) = rota_axis(3)
        valv(5) = rota_cent(1)
        valv(6) = rota_cent(2)
        valv(7) = rota_cent(3)

! ----- Read mesh affectation
        call getelem(mesh, keywordfact, iocc, ' ', listCell, nbCell)
        if (nbCell .eq. 0) then
            call nocart(map(1), 1, nbCmp(1))
        else
            call jeveuo(listCell, 'L', jvCell)
            call nocart(map(1), 3, nbCmp(1), mode='NUM', nma=nbCell, &
                        limanu=zi(jvCell))
        end if
!
        call jedetr(listCell)
    end do
!
99  continue
    call jedema()
end subroutine
