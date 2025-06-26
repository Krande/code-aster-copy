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
subroutine cavite(phenom, load, mesh, valeType, nbOcc)
!
    implicit none
!
#include "asterf_types.h"
#include "LoadTypes_type.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/char_crea_cart.h"
#include "asterfort/getelem.h"
#include "asterfort/getvc8.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
!
    character(len=16), intent(in) :: phenom
    character(len=8), intent(in) :: load, mesh
    character(len=4), intent(in) :: valeType
    integer(kind=8), intent(in) :: nbOcc
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation - Acoustic and mechanic
!
! Treatment of load VITE_FACE
!
! --------------------------------------------------------------------------------------------------
!
! In  phenom           : phenomenon (MECANIQUE/THERMIQUE/ACOUSTIQUE)
! In  load             : load
! In  mesh             : mesh
! In  valeType         : affected value type (real, complex or function)
! In  nbOcc            : number of factor keywords
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywFact = 'VITE_FACE'
    character(len=24), parameter :: listCell = '&&CAVITE.LIST_CELL'
    integer(kind=8) :: nbCell, jvCell
    integer(kind=8) :: iocc, nbRet, nbVal
    character(len=16) :: keyword
    complex(kind=8) :: speedCplx
    real(kind=8) :: speedReal, speedDirectionReal(3)
    character(len=8) :: speedFunc, speedDirectionFunc(3)
    integer(kind=8) :: jvValv
    character(len=19) :: map(LOAD_MAP_NBMAX)
    integer(kind=8) :: nbMap, nbCmp(LOAD_MAP_NBMAX)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Creation and initialization to zero of <CARTE>
    call char_crea_cart(phenom, keywFact, load, mesh, valeType, &
                        nbMap, map, nbCmp)
    ASSERT(nbMap .eq. 1)
    call jeveuo(map(1)//'.VALV', 'E', jvValv)

! - Loop on factor keyword
    do iocc = 1, nbOcc
! ----- Read mesh affectation
        call getelem(mesh, keywFact, iocc, 'A', listCell, nbCell)

! ----- Get type of speed: with direction or without it
        if (valeType .eq. 'REEL') then
            call getvr8(keywFact, 'VNOR', iocc=iocc, nbval=0, nbret=nbRet)
        elseif (valeType .eq. 'COMP') then
            call getvc8(keywFact, 'VNOR', iocc=iocc, nbval=0, nbret=nbRet)
            ASSERT(nbRet .eq. -1)
        elseif (valeType .eq. 'FONC') then
            call getvid(keywFact, 'VNOR', iocc=iocc, nbval=0, nbret=nbRet)
        else
            ASSERT(ASTER_FALSE)
        end if
        nbRet = abs(nbRet)
        ASSERT(nbRet .le. 1)
        if (nbRet .eq. 0) then
            keyword = 'VITE'
        else
            keyword = 'VNOR'
        end if
        if (valeType .eq. 'COMP') then
            ASSERT(keyword .eq. 'VNOR')
        end if

! ----- Get value of speed
        if (valeType .eq. 'REEL') then
            call getvr8(keywFact, keyword, iocc=iocc, scal=speedReal, nbret=nbRet)
        elseif (valeType .eq. 'COMP') then
            call getvc8(keywFact, keyword, iocc=iocc, scal=speedCplx, nbret=nbRet)
        elseif (valeType .eq. 'FONC') then
            call getvid(keywFact, keyword, iocc=iocc, scal=speedFunc, nbret=nbRet)
        end if

! ----- Get direction if necessary
        speedDirectionReal = 0.d0
        speedDirectionFunc = ' '
        if (keyword .eq. 'VITE') then
            if (valeType .eq. 'REEL') then
                call getvr8(keywFact, 'DIRECTION', &
                            iocc=iocc, nbval=0, nbret=nbRet)
                nbVal = abs(nbRet)
                call getvr8(keywFact, 'DIRECTION', &
                            iocc=iocc, nbval=nbVal, vect=speedDirectionReal, nbret=nbRet)
            elseif (valeType .eq. 'COMP') then
                speedDirectionReal = 0.d0
            elseif (valeType .eq. 'FONC') then
                call getvid(keywFact, 'DIRECTION', &
                            iocc=iocc, nbval=0, nbret=nbRet)
                nbVal = abs(nbRet)
                call getvid(keywFact, 'DIRECTION', &
                            iocc=iocc, nbval=nbVal, vect=speedDirectionFunc, nbret=nbRet)
            end if
        end if

! ----- Save values
        if (keyword .eq. 'VNOR') then
            if (valeType .eq. 'REEL') then
                zr(jvValv-1+1) = speedReal
                zr(jvValv-1+2) = -1.d0
            elseif (valeType .eq. 'COMP') then
                zc(jvValv-1+1) = speedCplx
                zc(jvValv-1+2) = (-1.d0, 0.d0)
            elseif (valeType .eq. 'FONC') then
                zk8(jvValv-1+1) = speedFunc
                zk8(jvValv-1+2) = '&FOZERO'
            end if
        else
            if (valeType .eq. 'REEL') then
                zr(jvValv-1+1) = speedReal
                zr(jvValv-1+2) = 1.d0
                zr(jvValv-1+3) = speedDirectionReal(1)
                zr(jvValv-1+4) = speedDirectionReal(2)
                zr(jvValv-1+5) = speedDirectionReal(3)
            elseif (valeType .eq. 'COMP') then
                zc(jvValv-1+1) = speedCplx
                zc(jvValv-1+2) = (1.d0, 0.d0)
                zc(jvValv-1+3) = speedDirectionReal(1)
                zc(jvValv-1+4) = speedDirectionReal(2)
                zc(jvValv-1+5) = speedDirectionReal(3)
            elseif (valeType .eq. 'FONC') then
                zk8(jvValv-1+1) = speedFunc
                zk8(jvValv-1+2) = ' '
                zk8(jvValv-1+3) = speedDirectionFunc(1)
                zk8(jvValv-1+4) = speedDirectionFunc(2)
                zk8(jvValv-1+5) = speedDirectionFunc(3)
            end if

        end if

! ----- Set parameter in field
        if (nbCell .ne. 0) then
            call jeveuo(listCell, 'L', jvCell)
            call nocart(map(1), 3, nbCmp(1), mode='NUM', nma=nbCell, limanu=zi(jvCell))
            call jedetr(listCell)
        end if

    end do
!
    call jedema()
end subroutine
