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
subroutine dmatdp(fami, materCodeJv, time, poum, ipg, &
                  ispg, anglNaut, dr_, di_)
!
    implicit none
!
#include "asterfort/get_elas_id.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/matrHookePlaneStrain.h"
#include "asterfort/separ_RI_elas_3D.h"
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: materCodeJv
    real(kind=8), intent(in) :: time
    character(len=*), intent(in) :: poum
    integer(kind=8), intent(in) :: ipg, ispg
    real(kind=8), intent(in) :: anglNaut(3)
    real(kind=8), optional, intent(out) :: dr_(4, 4), di_(4, 4)
!
! --------------------------------------------------------------------------------------------------
!
! Hooke matrix for iso-parametric elements
!
! Plane strain and axisymmetric
!
! --------------------------------------------------------------------------------------------------
!
! In  fami             : Gauss family for integration point rule
! In  materCodeJv      : coded material address
! In  time             : current time
! In  poum             : '-' or '+' for parameters evaluation (previous or current temperature)
! In  ipg              : current point gauss
! In  ispg             : current "sous-point" gauss
! In  anglNaut         : nautical angles for definition of basis for non-isotropic elasticity
! Out dr               : Hooke matrix (real part)
! Out di               : Hooke matrix (imaginary part, for viscoelasticity)
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: nur, nui, nu12r, nu13r, nu23r, nu12i, nu13i, nu23i
    real(kind=8) :: e1r, e2r, e3r, e1i, e2i, e3i, er, ei
    real(kind=8) :: g1r, g2r, g3r, g1i, g2i, g3i, gr, gi
    real(kind=8) :: di(4, 4), dr(4, 4), hr(6), hi(6)
    integer(kind=8) :: elasID
    character(len=32) :: elasKeyword
!
! --------------------------------------------------------------------------------------------------
!

! - Get type of elasticity (Isotropic/Orthotropic/Transverse isotropic)
    call get_elas_id(materCodeJv, elasID, elasKeyword)

! - Get elastic parameters
    call get_elas_para(fami, materCodeJv, poum, ipg, ispg, &
                       elasID, elasKeyword, &
                       time=time, &
                       e_=er, nu_=nur, g_=gr, &
                       e1_=e1r, e2_=e2r, e3_=e3r, &
                       nu12_=nu12r, nu13_=nu13r, nu23_=nu23r, &
                       g1_=g1r, g2_=g2r, g3_=g3r, &
                       ei_=ei, nui_=nui, gi_=gi, &
                       e1i_=e1i, e2i_=e2i, e3i_=e3i, &
                       nu12i_=nu12i, nu13i_=nu13i, nu23i_=nu23i, &
                       g1i_=g1i, g2i_=g2i, g3i_=g3i)

! - Compute Hooke matrix element
    call separ_RI_elas_3D(elasID, nur, gr, nui, gi, &
                          e1r, e2r, e3r, &
                          nu12r, nu13r, nu23r, &
                          e1i, e2i, e3i, &
                          nu12i, nu13i, nu23i, &
                          hr, hi)

! - Compute Hooke matrix
    if (present(di_)) then
        call matrHookePlaneStrain(elasID, anglNaut, &
                                  hi, gi, g1i, &
                                  di)
        di_ = di
    end if
    if (present(dr_)) then
        call matrHookePlaneStrain(elasID, anglNaut, &
                                  hr, gr, g1r, &
                                  dr)
        dr_ = dr
    end if
!
end subroutine
