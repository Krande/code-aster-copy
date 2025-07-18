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

subroutine lcmfra(vp, itemax, precvg, chi, iret)
    implicit none
#include "asterc/r8prem.h"
    integer(kind=8), intent(in)      :: itemax
    real(kind=8), intent(in) :: vp(3), precvg
    integer(kind=8), intent(out)     :: iret
    real(kind=8), intent(out):: chi
! --------------------------------------------------------------------------------------------------
!  CRITERE ENDO_FISS_EXP: CALCUL DU RAYON SEUIL
! --------------------------------------------------------------------------------------------------
! IN  VP      CONTRAINTES PRINCIPALES RANGEES DANS L'ORDRE DECROISSANT
! IN  ITEMAX  NOMBRE MAX D'ITERATIONS POUR LA RESOLUTION SCALAIRE
! IN  PRECVG  PRECISION DE LA RESOLUTION SCALAIRE : CHI*DCHI < PRECVG
! OUT CHI     SCALING TEL QUE VP/CHI SOIT SUR LA FRONTIERE DU CRITERE
! OUT ITET    CODE RETOUR: 0=OK, 1=NON CONVERGENCE
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: iter
    real(kind=8) :: nors, n(3), tra, nora, ymax, y, e(3), nore, f, df, gamma, prectr
    real(kind=8) :: precy3, yt, ft
! --------------------------------------------------------------------------------------------------
    real(kind=8) :: tau, sig0, beta
    common/lcmmf/tau, sig0, beta
! --------------------------------------------------------------------------------------------------
!
!  INITIALISATION
    iret = 0
    prectr = 1.d-3*tau
    gamma = 3*beta**2-1.d0/3.d0
!
!
!  NORMALISATION
    nors = sqrt((vp(1)**2+vp(2)**2+vp(3)**2))
    if (nors .lt. sig0*r8prem()) then
        chi = 0
        goto 999
    end if
!
    n(1) = vp(1)/nors
    n(2) = vp(2)/nors
    n(3) = vp(3)/nors
!
!
!  COEFFICIENT LINEAIRE
    tra = n(1)+n(2)+n(3)
    nora = sqrt(1+gamma*tra**2)
!
!
!  BORNE MAX
    ymax = tau/nora
    y = ymax
    if (n(1) .gt. 0) y = min(y, log(tau)/n(1))
!
!
!  RESOLUTION PAR UNE METHODE DE NEWTON

!  PRECISION REQUISE
    precy3 = (sig0/nors)**2*precvg
    do iter = 1, itemax
        e = exp(2*y*n)
        nore = sqrt(sum(e))
!
        f = nora*y+nore-tau
        df = nora+dot_product(n, e)/nore
!
        if (abs(f) .le. prectr) then
            yt = y-sign(1.d0, f)*precy3*y**3
            if (yt <= 0 .or. yt >= ymax) goto 20
            ft = nora*yt+sqrt(sum(exp(2*yt*n)))-tau
            if (ft*f .le. 0) goto 20
        end if

        y = y-f/df
    end do
!
!  ECHEC DANS LA RESOLUTION
    iret = 1
    goto 999
!
!  SUCCES DE LA RESOLUTION
20  continue
    chi = nors/(y*sig0)
!
!
999 continue
end subroutine
