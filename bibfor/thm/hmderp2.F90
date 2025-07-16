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
! aslint: disable=W1504
!
subroutine hmderp2(ds_thm, t, pvp, &
                   rho11, rho12, h11, h12, &
                   dp11p1, dp11p2, dp11t, &
                   dp12p1, dp12p2, dp12t, &
                   dp21p1, dp21p2, dp21t, &
                   dp1pp1, dp2pp1, dtpp1, &
                   dp1pp2, dp2pp2, dtpp2, &
                   dp1pt, dp2pt, dtpt)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
!
    type(THM_DS), intent(in) :: ds_thm
    real(kind=8), intent(in) :: t, pvp
    real(kind=8), intent(in) :: rho11, rho12, h11, h12
    real(kind=8), intent(out) :: dp11p1, dp11p2, dp11t
    real(kind=8), intent(out) :: dp12p1, dp12p2, dp12t
    real(kind=8), intent(out) :: dp21p1, dp21p2, dp21t
    real(kind=8), intent(out) :: dp1pp1, dp2pp1, dtpp1
    real(kind=8), intent(out) :: dp1pp2, dp2pp2, dtpp2
    real(kind=8), intent(out) :: dp1pt, dp2pt, dtpt
!
! -------------------------------------------------------------------------------------------------
!
! THM
!
! Compute some derivatives for LIQU_VAPE_GAZ
!
! -------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! In  t                : temperature - At end of current step
! In  pvp              : steam pressure
! In  pad              : dissolved air pressure
! In  rho11            : current volumic mass of liquid
! In  rho12            : current volumic mass of steam
! In  h11              : enthalpy of liquid
! In  h12              : enthalpy of steam
! Out dp11p1           : derivative of P11 by P1
! Out dp11p2           : derivative of P11 by P2
! Out dp11t            : derivative of P11 by T
! Out dp12p1           : derivative of P12 by P1
! Out dp12p2           : derivative of P12 by P2
! Out dp12t            : derivative of P12 by T
! Out dp21p1           : derivative of P21 by P1
! Out dp21p2           : derivative of P21 by P2
! Out dp21t            : derivative of P21 by T
! Out dp1pp1           : second derivative of P1 by P1²
! Out dp2pp1           : second derivative of P2 by P1²
! Out dtpp1            : second derivative of T by P1²
! Out dp1pp2           : second derivative of P1 by P2²
! Out dp2pp2           : second derivative of P2 by P2²
! Out dtpp2            : second derivative of T by P2²
! Out dp1pt            : derivative of P1 by T
! Out dp2pt            : derivative of P2 by T
! Out dtpt             : derivative of T by T
!
! -------------------------------------------------------------------------------------------------
!
    real(kind=8) :: cliq, alpliq
    real(kind=8) :: a1, a3, a4, l
    real(kind=8), parameter :: zero = 0.d0
!
! -------------------------------------------------------------------------------------------------
!
    dp11p1 = zero
    dp11p2 = zero
    dp11t = zero
    dp12p1 = zero
    dp12p2 = zero
    dp12t = zero
    dp21p1 = zero
    dp21p2 = zero
    dp21t = zero
    dp1pp1 = zero
    dp2pp1 = zero
    dtpp1 = zero
    dp1pp2 = zero
    dp2pp2 = zero
    dtpp2 = zero
    dp1pt = zero
    dp2pt = zero
    dtpt = zero
!
! - Get parameters

    cliq = ds_thm%ds_material%liquid%unsurk
    alpliq = ds_thm%ds_material%liquid%alpha
!
! - Compute
!
    dp11p1 = -1.
    dp11p2 = 1.
    dp12p1 = -rho12/rho11
    dp12p2 = rho12/rho11
    dp21p1 = -dp12p1
    dp21p2 = 1.d0-dp12p2
    if (ds_thm%ds_elem%l_dof_ther) then
        dp11t = 0.d0
        dp12t = rho12*(h12-h11)/t
        dp21t = -dp12t
        l = (h12-h11)
    end if
    a1 = -rho11/rho12
    if (ds_thm%ds_elem%l_dof_ther) then
        a3 = dp12t/pvp-1/t-3*alpliq
        a4 = -dp12t/pvp+1/t-3*alpliq
    else
        a3 = dp12t/pvp-1/t
        a4 = -dp12t/pvp+1/t
    end if
    dp1pp2 = (cliq*dp11p1-1/pvp*dp12p1)/a1
    dp2pp2 = (cliq*dp11p2-1/pvp*dp12p2)/a1
    dp1pp1 = -rho11/rho12/a1/a1*(dp12p1/pvp-cliq*dp11p1)
    dp2pp1 = -rho11/rho12/a1/a1*(dp12p2/pvp-cliq*dp11p2)
    if (ds_thm%ds_elem%l_dof_ther) then
        dtpp2 = 1./a1/a1*(-a4*rho11/rho12)
        dtpp1 = -1.d0/a1/a1*(-rho11/rho12*a4)
        dp1pt = -1.d0/t/a1/a1*(a1*(1.d0-dp11p1*(1.d0+l*rho11*cliq)) &
                               -l*rho11*rho11/rho12*(cliq*dp11p1-dp12p1/pvp))
        dp2pt = -1.d0/t/a1/a1*(a1*(1.d0-dp11p2*(1.d0+l*rho11*cliq)) &
                               -l*rho11*rho11/rho12*(cliq*dp11p2-dp12p2/pvp))
        dtpt = -1.d0/t/a1*l*(rho11*dp11t*cliq-3.d0*alpliq*rho11)- &
               1.d0/t/t/a1/a1*(rho11/rho12+t*(rho11/rho12*a4))*(l*rho11)

    end if
!
end subroutine
