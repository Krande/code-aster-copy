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

subroutine te0083(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
!
    character(len=16) :: option, nomte
!     OPTION : EFEQ_ELNO
!     IN   K16   OPTION : NOM DE L'OPTION A CALCULER
!     IN   K16   NOMTE  : NOM DU TYPE_ELEMENT
!
!     POST_RCCM/MOMENT_EQUIVALENT
!     4 composantes en sortie :
!     MT, MFY, MFZ et MEQ = sqrt(MT**2+MFY**2+MFZ**2)
!     ------------------------------------------------------------------
    integer(kind=8) :: jin, jout, i, ibid, itab(7), lgcata
    integer(kind=8) :: nbpoin, nbcompin, nbcompout
    real(kind=8) :: mt, mfy, mfz, meq
!     ------------------------------------------------------------------
!
    ASSERT(option .eq. 'EFEQ_ELNO')
    call jevech('PEFFONR', 'L', jin)
    call tecach('OOO', 'PEFFONR', 'L', ibid, nval=7, &
                itab=itab)
    call jevech('PEFFOENR', 'E', jout)

    jin = itab(1)
    nbpoin = itab(3)
    lgcata = itab(2)
    nbcompin = lgcata/nbpoin
    ASSERT(nbcompin .eq. 6 .or. nbcompin .eq. 7)
    nbcompout = 4
    do i = 1, nbpoin
        mt = zr(jin+nbcompin*(i-1)-1+4)
        mfy = zr(jin+nbcompin*(i-1)-1+5)
        mfz = zr(jin+nbcompin*(i-1)-1+6)
        meq = sqrt(mt*mt+mfy*mfy+mfz*mfz)
        zr(jout+nbcompout*(i-1)-1+1) = mt
        zr(jout+nbcompout*(i-1)-1+2) = mfy
        zr(jout+nbcompout*(i-1)-1+3) = mfz
        zr(jout+nbcompout*(i-1)-1+4) = meq
    end do
!
end subroutine
