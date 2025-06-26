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
subroutine inmat4(elrefa, nno, nnos, npg, nofpg, mgano)
!
    implicit none
!
#include "MeshTypes_type.h"
#include "asterfort/assert.h"
#include "asterfort/inmat5.h"
#include "asterfort/inmat6.h"
!
    character(len=8) :: elrefa, nofpg
    integer(kind=8) :: nno, nnos, npg
    real(kind=8) :: mgano(*)

! ======================================================================
! BUT : CALCULER LA MATRICE DE PASSAGE GAUSS -> NOEUDS
!       POUR UNE FAMILLE D'UN ELREFA
! ======================================================================
!
    integer(kind=8) :: kpg, kno, knos, k
    real(kind=8) :: mganos(MT_NBPGMX, MT_NNOMAX), mgano2(MT_NBPGMX, MT_NNOMAX)
    ASSERT(npg .le. MT_NBPGMX)
    ASSERT(nno .le. MT_NNOMAX)
    ASSERT(nnos .le. MT_NNOMAX)
!
!
!     -- MISES A ZERO :
!     ----------------------------------------------------------
    do kpg = 1, npg
        do kno = 1, nno
            mgano2(kpg, kno) = 0.d0
        end do
        do knos = 1, nnos
            mganos(kpg, knos) = 0.d0
        end do
    end do
    do k = 1, 2+npg*nno
        mgano(k) = 0.d0
    end do
!
!
!     -- ON TRAITE LE CAS GENERIQUE NPG=1  (INCLUT NOFPG='FPG1')
!     ----------------------------------------------------------
    if (npg .eq. 1) then
        do kno = 1, nno
            mgano2(1, kno) = 1.d0
        end do
        goto 80
    end if
!
!
!     -- ON TRAITE LE CAS GENERIQUE NOFPG='NOEU'
!     -------------------------------------------------
    if (nofpg .eq. 'NOEU') then
        ASSERT(nno .eq. npg)
        do k = 1, nno
            mgano2(k, k) = 1.d0
        end do
        goto 80
    end if
!
!
!     -- ON TRAITE LE CAS GENERIQUE NOFPG='NOEU_S'
!     -------------------------------------------------
    if (nofpg .eq. 'NOEU_S') then
        ASSERT(nnos .eq. npg)
        do k = 1, nnos
            mganos(k, k) = 1.d0
        end do
        call inmat5(elrefa, nno, nnos, npg, mganos, &
                    mgano2)
        goto 80
    end if
!
!
!     -- AUTRES CAS : GAUSS -> SOMMETS -> NOEUDS
!     -------------------------------------------
    call inmat6(elrefa, nofpg, mganos)
    call inmat5(elrefa, nno, nnos, npg, mganos, &
                mgano2)
!
80  continue
    mgano(1) = nno
    mgano(2) = npg
    do kpg = 1, npg
        do kno = 1, nno
            mgano(2+(kno-1)*npg+kpg) = mgano2(kpg, kno)
        end do
    end do
!
end subroutine
