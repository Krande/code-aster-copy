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
subroutine caelec(load, mesh, nbOcc)
!
    implicit none
!
#include "asterf_types.h"
#include "LoadTypes_type.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/char_crea_cart.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/getelem.h"
!
    character(len=8), intent(in) :: load, mesh
    integer(kind=8), intent(in) :: nbOcc
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load FORCE_ELEC
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  mesh             : mesh
! In  nbOcc            : number of factor keywords
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordfact = 'FORCE_ELEC'
    character(len=4), parameter :: valeType = 'REEL'
    character(len=24), parameter :: listCell = '&&CAELEC.LIST_ELEM'
    integer(kind=8) :: jvCell, nbCell
    integer(kind=8) :: nbRet, iocc
    real(kind=8) :: p1(3), p2(3), zcod, d
    character(len=8) :: code
    real(kind=8), pointer :: valv(:) => null()
    character(len=19) :: map(LOAD_MAP_NBMAX)
    integer(kind=8) :: nbMap, nbCmp(LOAD_MAP_NBMAX)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
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
    do iocc = 1, nbOcc

! ----- Get values
        call getvtx(keywordfact, 'POSITION', iocc=iocc, scal=code, nbret=nbRet)
        p1 = 0.d0
        p2 = 0.d0
        if (nbRet .eq. 0) then
            zcod = 10.d0
            call getvr8(keywordfact, 'FX', iocc=iocc, scal=p1(1), nbret=nbRet)
            call getvr8(keywordfact, 'FY', iocc=iocc, scal=p1(2), nbret=nbRet)
            call getvr8(keywordfact, 'FZ', iocc=iocc, scal=p1(3), nbret=nbRet)
        else
            if (code .eq. 'PARA') then
                call getvr8(keywordfact, 'DIST', iocc=iocc, scal=d, nbret=nbRet)
                if (nbRet .ne. 0) then
                    zcod = 12.d0
                    p1(1) = d
                    p1(2) = 0.d0
                    p1(3) = 0.d0
                    call getvr8(keywordfact, 'POINT2', iocc=iocc, nbval=3, vect=p2, nbret=nbRet)
                else
                    zcod = 11.d0
                    call getvr8(keywordfact, 'TRANS', iocc=iocc, nbval=3, vect=p1, nbret=nbRet)
                    p2(1) = 0.d0
                    p2(2) = 0.d0
                    p2(3) = 0.d0
                end if
            else if (code .eq. 'INFI') then
                zcod = 2.d0
                call getvr8(keywordfact, 'POINT1', iocc=iocc, nbval=3, vect=p1, nbret=nbRet)
                call getvr8(keywordfact, 'POINT2', iocc=iocc, nbval=3, vect=p2, nbret=nbRet)
            else if (code .eq. 'FINI') then
                zcod = 3.d0
                call getvr8(keywordfact, 'POINT1', iocc=iocc, nbval=3, vect=p1, nbret=nbRet)
                call getvr8(keywordfact, 'POINT2', iocc=iocc, nbval=3, vect=p2, nbret=nbRet)
            end if
        end if

! ----- Read mesh affectation
        call getelem(mesh, keywordfact, iocc, 'F', listCell, nbCell)
        call jeveuo(listCell, 'L', jvCell)

! ----- Set values in map
        valv(1) = p1(1)
        valv(2) = p1(2)
        valv(3) = p1(3)
        valv(4) = p2(1)
        valv(5) = p2(2)
        valv(6) = p2(3)
        valv(7) = zcod
        call nocart(map(1), 3, nbCmp(1), mode='NUM', nma=nbCell, &
                    limanu=zi(jvCell))
        call jedetr(listCell)
    end do
!
    call jedema()
end subroutine
