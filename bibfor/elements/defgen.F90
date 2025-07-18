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
subroutine defgen(testl1, testl2, nno, r, x3, &
                  sina, cosa, cour, vf, dfds, &
                  depl, eps, epsx3)
    implicit none
#include "asterf_types.h"
!
    integer(kind=8) :: nno
    aster_logical :: testl1, testl2
    real(kind=8) :: r, x3, sina, cosa, cour, vf(*), dfds(*), depl(*), eps(*)
    real(kind=8) :: epsx3
    real(kind=8) :: uxl(3), uyl(3), betasl(3)
    real(kind=8) :: ess, kss, ett, ktt, gs
!
!     CALCUL DES DEFORMATIONS GENERALISEES : ESS , ETT , KSS , KTT , GS
!
!     PARTITION DU DEPL EN UX, UY ET BETAS
!
!     DO 10 INO=1,NNO
!-----------------------------------------------------------------------
    integer(kind=8) :: i
    real(kind=8) :: betas, dbtds, duxds, duyds, rhos, rhot, ux
    real(kind=8) :: uy
!-----------------------------------------------------------------------
    uxl(1) = depl(1)
    uxl(2) = depl(4)
    uxl(3) = depl(7)
!
    uyl(1) = depl(2)
    uyl(2) = depl(5)
    uyl(3) = depl(8)
!
    betasl(1) = depl(3)
    betasl(2) = depl(6)
    betasl(3) = depl(9)
!10   CONTINUE
!
    ux = 0.d0
    uy = 0.d0
    betas = 0.d0
!
    duxds = 0.d0
    duyds = 0.d0
    dbtds = 0.d0
    do i = 1, nno
        ux = ux+vf(i)*uxl(i)
        uy = uy+vf(i)*uyl(i)
        betas = betas+vf(i)*betasl(i)
!
        duxds = duxds+dfds(i)*uxl(i)
        duyds = duyds+dfds(i)*uyl(i)
        dbtds = dbtds+dfds(i)*betasl(i)
    end do
!
!     ESS  ,  KSS  ,  ETT  ,  KTT  ,  GS
!
    ess = duyds*cosa-duxds*sina
    kss = dbtds
    ett = ux/r
    ktt = -sina/r*betas
    gs = betas+duxds*cosa+duyds*sina
!
    if (testl1) then
        rhos = 1.d0
    else
        rhos = 1.d0+x3*cour
    end if
    if (testl2) then
        rhot = 1.d0
    else
        rhot = 1.d0+x3*cosa/r
    end if
!
    eps(1) = (ess+x3*kss)/rhos
    eps(2) = (ett+x3*ktt)/rhot
    eps(3) = 0.d0
    eps(4) = 0.d0
!
    epsx3 = 0.5d0/rhos*gs
!
end subroutine
