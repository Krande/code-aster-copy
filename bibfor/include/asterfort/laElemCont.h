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
#include "asterf_types.h"
#include "contact_module.h"
!
interface
    subroutine laElemCont(parameters, geom, coor_qp_sl, hF, &
                          gap, gamma_c, projRmVal, l_cont_qp, &
                          gamma_f, l_fric_qp, &
                         norm_slav, dGap, projBsVal, dv, dvT, projBsVal2, lagr_c, lagr_f, lagr_f3, &
                          dThres_du, dThres_dl, lagr_v, thres, jump_v, dNs, &
                          tau_slav, lagr_g, mu_g, dts_ns, speed, dfunc_ma, dZetaM, &
                          mu_c, mu_f, mu_f3, metricTens, d2Gap)
!
        use contact_type
!
        type(ContactParameters), intent(in) :: parameters
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: coor_qp_sl(2), hF
        real(kind=8), intent(out) :: gap, gamma_c, projRmVal
        real(kind=8), intent(out) :: gamma_f
        aster_logical, intent(out) :: l_cont_qp, l_fric_qp
        real(kind=8), intent(out), optional :: dts_ns(MAX_LAGA_DOFS, 2), dGap(MAX_LAGA_DOFS), thres
        real(kind=8), intent(out), optional :: speed(3), dvT(MAX_LAGA_DOFS, 2)
        real(kind=8), intent(out), optional :: dfunc_ma(3, MAX_LAGA_DOFS, 2)
        real(kind=8), intent(out), optional :: mu_f(MAX_LAGA_DOFS, 2), mu_c(MAX_LAGA_DOFS)
        real(kind=8), intent(out), optional :: tau_slav(3, 2), projBsVal2(2), lagr_c, lagr_f(2)
        real(kind=8), intent(out), optional :: lagr_g(3), mu_g(MAX_LAGA_DOFS, 3), metricTens(2, 2)
        real(kind=8), intent(out), optional :: jump_v(MAX_LAGA_DOFS, 3), dZetaM(MAX_LAGA_DOFS, 2)
        real(kind=8), intent(out), optional :: projBsVal(3), lagr_v(3), dNs(MAX_LAGA_DOFS, 3)
        real(kind=8), intent(out), optional :: dv(MAX_LAGA_DOFS, 3), norm_slav(3)
        real(kind=8), intent(out), optional :: d2Gap(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
        real(kind=8), intent(out), optional :: dThres_du(MAX_LAGA_DOFS), dThres_dl(MAX_LAGA_DOFS)
        real(kind=8), intent(out), optional :: lagr_f3(3), mu_f3(MAX_LAGA_DOFS, 3)
!
    end subroutine laElemCont
end interface
