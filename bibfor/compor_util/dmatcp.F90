! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

subroutine dmatcp(fami, mater, time, poum, ipg,&
                  ispg, repere, dr_, di_)
!
implicit none
!
#include "asterfort/assert.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/matrHookePlaneStress.h"
!
!
    character(len=*), intent(in) :: fami
    integer, intent(in) :: mater
    real(kind=8), intent(in) :: time
    character(len=*), intent(in) :: poum
    integer, intent(in) :: ipg
    integer, intent(in) :: ispg
    real(kind=8), intent(in) :: repere(7)
    real(kind=8), optional, intent(out) :: dr_(4, 4)
    real(kind=8), optional, intent(out) :: di_(4, 4)
!
! --------------------------------------------------------------------------------------------------
!
! Hooke matrix for iso-parametric elements
!
! Plane stress
!
! --------------------------------------------------------------------------------------------------
!
! In  fami   : Gauss family for integration point rule
! In  mater  : material parameters
! In  time   : current time
! In  poum   : '-' or '+' for parameters evaluation (previous or current temperature)
! In  ipg    : current point gauss
! In  ispg   : current "sous-point" gauss
! In  repere : local basis for orthotropic elasticity
! Out dr     : real Hooke matrix
! Out di     : imaginary Hooke matrix
!
! --------------------------------------------------------------------------------------------------
!
    integer :: elas_id
    real(kind=8) :: nur, nui, nu12r, nu13r, nu23r, nu12i, nu13i, nu23i
    real(kind=8) :: e1r, e2r, e3r, e1i, e2i, e3i, er, ei
    real(kind=8) :: g1r, g2r, g3r, g1i, g2i, g3i, gr, gi
    character(len=16) :: elas_keyword
    real(kind=8) :: di(4, 4), dr(4,4)
!
! --------------------------------------------------------------------------------------------------
!
!
! - Get type of elasticity (Isotropic/Orthotropic/Transverse isotropic)
!
    call get_elas_id(mater, elas_id, elas_keyword)
!
! - Get elastic parameters
!
    call get_elas_para(fami, mater    , poum, ipg, ispg, &
                       elas_id  , elas_keyword,&
                       time = time,&
                       e_ = er    , nu_ = nur  , g_ = gr,&
                       e1_ = e1r    , e2_ = e2r    , e3_ = e3r,&
                       nu12_ = nu12r, nu13_ = nu13r, nu23_ = nu23r,&
                       g1_ = g1r    , g2_ = g2r    , g3_ = g3r,&
                       ei_ = ei    , nui_ = nui  , gi_ = gi,&
                       e1i_ = e1i    , e2i_ = e2i    , e3i_ = e3i,&
                       nu12i_ = nu12i, nu13i_ = nu13i, nu23i_ = nu23i,&
                       g1i_ = g1i    , g2i_ = g2i    , g3i_ = g3i)
!
! - Compute Hooke matrix
!
    if (present(di_)) then
        call matrHookePlaneStress(elas_id, repere,&
                                  ei , nui,&
                                  e1i, e2i, nu12i, g1i,&
                                  di)
        di_ = di
    endif
    if (present(dr_)) then
        call matrHookePlaneStress(elas_id, repere,&
                                  er , nur,&
                                  e1r, e2r, nu12r, g1r,&
                                  dr)
        dr_ = dr
    endif
!
end subroutine
