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

subroutine xellin(elref1, nno1, elref2, nno2)
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/iselli.h"
    character(len=8) :: elref1, elref2
    integer(kind=8) :: nno1, nno2
!
!                      RETOURNE LE TYPE DE L'ELEMENT "LINEARISE"
!                      ET LE NOMBRE DE NOEUDS DE CHAQUE ELEMENT
!
!     ENTREE
!       ELREF1  : TYPE DE L'ELEMENT PARENT
!       NNO1    : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
!
!     SORTIE
!       ELREF2  : TYPE DE L'ELEMENT LINEAIRE A L'ELEMENT PARENT

!       NNO2    : NOMBRE DE NOEUDS DE L'ELEMENT LINEAIRE
!......................................................................
!
    call jemarq()
!
    if (iselli(elref1)) goto 999
!
    if (elref1 .eq. 'QU8') then
        ASSERT(nno1 .eq. 8)
        elref2 = 'QU4'
        nno2 = 4
    else if (elref1 .eq. 'TR6') then
        ASSERT(nno1 .eq. 6)
        elref2 = 'TR3'
        nno2 = 3
    else if (elref1 .eq. 'SE3') then
        ASSERT(nno1 .eq. 3)
        elref2 = 'SE2'
        nno2 = 2
    else if (elref1 .eq. 'H20') then
        ASSERT(nno1 .eq. 20)
        elref2 = 'HE8'
        nno2 = 8
    else if (elref1 .eq. 'P15') then
        ASSERT(nno1 .eq. 15)
        elref2 = 'PE6'
        nno2 = 6
    else if (elref1 .eq. 'P13') then
        ASSERT(nno1 .eq. 13)
        elref2 = 'PY5'
        nno2 = 5
    else if (elref1 .eq. 'T10') then
        ASSERT(nno1 .eq. 10)
        elref2 = 'TE4'
        nno2 = 4
    else
        ASSERT(.false.)
    end if
!
999 continue
!
    call jedema()
!
end subroutine
