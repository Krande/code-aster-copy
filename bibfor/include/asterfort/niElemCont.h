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
    subroutine niElemCont(parameters, geom, nits, coor_qp_sl, hF, &
                    stress_nn, gap, gamma_c, projRmVal, l_cont_qp,&
                    stress_t, vT, gamma_f, projBsVal, l_fric_qp, &
                    dGap, d2Gap, jump_t, dStress_nn)
        use contact_type
        use contact_nitsche_module
        type(ContactParameters), intent(in) :: parameters
        type(ContactGeom), intent(in) :: geom
        type(ContactNitsche), intent(in) :: nits
        real(kind=8), intent(in) :: coor_qp_sl(2), hF
        real(kind=8), intent(out) :: stress_nn, gap, gamma_c, projRmVal
        real(kind=8), intent(out) :: stress_t(2), vT(2), gamma_f, projBsVal(2)
        aster_logical, intent(out) :: l_cont_qp, l_fric_qp
        real(kind=8), intent(out), optional :: dGap(MAX_LAGA_DOFS)
        real(kind=8), intent(out), optional :: d2Gap(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
        real(kind=8), intent(out), optional :: jump_t(MAX_LAGA_DOFS,3)
        real(kind=8), intent(out), optional :: dStress_nn(MAX_NITS_DOFS)
    end subroutine niElemCont
end interface
