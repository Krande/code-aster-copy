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

subroutine angvxz(gx, gn, angl)
    implicit none
!
    real(kind=8) :: gx(3), gn(3), angl(*)
!
! --------------------------------------------------------------------------------------------------
!
!   Calcul des 3 angles nautiques à partir du vecteur Gx et d'un vecteur Gn dont la projection
!   normale sur le plan normal à Gx donne le vecteur Gz
!
! --------------------------------------------------------------------------------------------------
!
!   IN  :   gx , gn
!   OUT :   angl
!
! --------------------------------------------------------------------------------------------------
!
#include "asterc/r8miem.h"
#include "asterc/r8pi.h"
#include "asterfort/angvx.h"
#include "asterfort/matrot.h"
#include "asterfort/pmavec.h"
!
    real(kind=8) :: mro(3, 3), gz(3)
    real(kind=8) :: alpha, beta, tst
!-----------------------------------------------------------------------
    tst = r8miem()
!
    call angvx(gx, alpha, beta)
    angl(1) = alpha
    angl(2) = beta
    angl(3) = 0.d0
    call matrot(angl, mro)
    call pmavec('ZERO', 3, mro, gn, gz)
    if ((abs(gz(3)) .le. tst) .and. (abs(gz(2)) .le. tst)) then
        angl(3) = r8pi()/2.0d0
    else
        angl(3) = atan2(-gz(2), gz(3))
    end if
!
end subroutine
