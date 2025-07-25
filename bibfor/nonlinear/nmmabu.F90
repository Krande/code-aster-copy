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
subroutine nmmabu(ndim, nno, axi, grand, dfdi, &
                  b)
!
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/r8inir.h"
#include "asterfort/utmess.h"
    aster_logical :: grand, axi
    integer(kind=8) :: ndim, nno
    real(kind=8) :: dfdi(nno, ndim), b(6, 3, nno)
!
! ----------------------------------------------------------------------
!                     CALCUL DE LA MATRICE B :  DEPS = B.DU
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  NNO     : NOMBRE DE NOEUDS
! IN  AXI     : .TRUE. SI AXISYMETRIQUE
! IN  GRAND   : .TRUE. SI GRANDES DEFORMATIONS
! IN  DFDI    : DERIVEE DES FONCTIONS DE FORME (POINT DE GAUSS COURANT)
! OUT B       : MATRICE B : B(6,3,NNP)
! ----------------------------------------------------------------------
!
    integer(kind=8) :: n
    real(kind=8) :: r2
! ----------------------------------------------------------------------
!
!
    if (grand) then
        call utmess('F', 'ALGORITH7_76')
    end if
    if (axi) then
        call utmess('F', 'ALGORITH7_76')
    end if
!
    call r8inir(18*nno, 0.d0, b, 1)
    r2 = sqrt(2.d0)/2.d0
!
    ASSERT((ndim .eq. 2) .or. (ndim .eq. 3))
!
    if (ndim .eq. 2) then
        do n = 1, nno
            b(1, 1, n) = dfdi(n, 1)
            b(2, 2, n) = dfdi(n, 2)
            b(4, 1, n) = r2*dfdi(n, 2)
            b(4, 2, n) = r2*dfdi(n, 1)
        end do
!
    else if (ndim .eq. 3) then
        do n = 1, nno
            b(1, 1, n) = dfdi(n, 1)
            b(2, 2, n) = dfdi(n, 2)
            b(3, 3, n) = dfdi(n, 3)
            b(4, 1, n) = r2*dfdi(n, 2)
            b(4, 2, n) = r2*dfdi(n, 1)
            b(5, 1, n) = r2*dfdi(n, 3)
            b(5, 3, n) = r2*dfdi(n, 1)
            b(6, 2, n) = r2*dfdi(n, 3)
            b(6, 3, n) = r2*dfdi(n, 2)
        end do
!
    end if
!
end subroutine
