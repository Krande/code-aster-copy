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

subroutine lcmfdr(sig, vp, chi, precvg, dchids)
    implicit none
#include "asterfort/lcexpo.h"
    real(kind=8), intent(in) :: sig(6)
    real(kind=8), intent(in) :: vp(3)
    real(kind=8), intent(in) :: chi
    real(kind=8), intent(in) :: precvg
    real(kind=8), intent(out):: dchids(6)
! --------------------------------------------------------------------------------------------------
!  CRITERE ENDO_FISS_EXP: DERIVEE DU FACTEUR D'ECHELLE SEUIL (CHI)
! sig(6)    contraintes
! vp(3)     contraintes principales
! chi       scaling tel que sig/chi soit sur la frontiere du critere
! precvg    précision du calcul (pour l'exponentiel matriciel)
! dchids(6) dérivée de chi par rapport a sig (1:6)
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter:: mexp = 6, factom = 720
    real(kind=8), parameter, dimension(6):: kr = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    real(kind=8), parameter :: srac2 = sqrt(0.5d0)
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: nexp
    real(kind=8) :: sbnor, sbtr, gamma, sbtnor, u(6), preexp, vpmax, targ
    real(kind=8) :: sb(6), e(6), e2(6), enor, h
! --------------------------------------------------------------------------------------------------
    real(kind=8) :: tau, sig0, beta
    common/lcmmf/tau, sig0, beta
! --------------------------------------------------------------------------------------------------
!
!  INITIALISATION
    preexp = precvg*1.d-1
    gamma = 3*beta**2-1.d0/3.d0
    sb = sig/sig0
    sbnor = sqrt(dot_product(sb, sb))
    u = sig/chi/sig0
!
!
!  CALCUL DE N TILDA
    sbtr = sb(1)+sb(2)+sb(3)
    sbtnor = sqrt(gamma*sbtr**2+sbnor**2)
!
!

!  CALCUL DE EXP(U) ET SA NORME
    vpmax = maxval(abs(vp/chi/sig0))
    targ = (factom*preexp)**(1.d0/mexp)
    nexp = merge(0, int(1+log(vpmax/targ)/log(2.d0)/(1-1.d0/mexp)), vpmax .lt. targ)
    call lcexpo(u, e, mexp, nexp)
    enor = sqrt(dot_product(e, e))

!  Exponentielle negligeable ou non
    if (enor*sbnor .le. sbtnor*preexp*1.d-3) then

        dchids = chi/sbtnor**2/sig0*(sb+gamma*sbtr*kr)

    else

!      CALCUL DE EXP(2.U)
        e2(1) = e(1)*e(1)+0.5d0*e(4)*e(4)+0.5d0*e(5)*e(5)
        e2(2) = 0.5d0*e(4)*e(4)+e(2)*e(2)+0.5d0*e(6)*e(6)
        e2(3) = 0.5d0*e(5)*e(5)+0.5d0*e(6)*e(6)+e(3)*e(3)
        e2(4) = e(1)*e(4)+e(4)*e(2)+srac2*e(5)*e(6)
        e2(5) = e(1)*e(5)+srac2*e(4)*e(6)+e(5)*e(3)
        e2(6) = srac2*e(4)*e(5)+e(2)*e(6)+e(6)*e(3)
!
!      CALCUL DE H
        h = sbtnor+sum(e2*sb)/enor
!
!      CALCUL DU TENSEUR TANGENT
        dchids = chi/h/sig0*((sb+gamma*sbtr*kr)/sbtnor+e2/enor)

    end if
!
end subroutine
