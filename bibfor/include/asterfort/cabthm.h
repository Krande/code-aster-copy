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
#include "asterf_types.h"
!
interface
    subroutine cabthm(ds_thm   , l_axi    , ndim   ,&
                      nddls    , nddlm    ,&
                      nddl_meca, nddl_p1  , nddl_p2, nddl_2nd, &
                      nno      , nnos     , &
                      dimuel   , dimdef   , kpi    ,&
                      addeme   , addete   , addep1 , addep2, adde2nd, &
                      elem_coor,&
                      jv_poids , jv_poids2,&
                      jv_func  , jv_func2 ,&
                      jv_dfunc , jv_dfunc2,&
                      dfdi     , dfdi2    ,&
                      poids    , poids2   ,&
                      b        )
        use THM_type
        type(THM_DS), intent(in) :: ds_thm
        aster_logical, intent(in) :: l_axi
        integer(kind=8), intent(in) :: ndim, nddls, nddlm
        integer(kind=8), intent(in) :: nddl_meca, nddl_p1, nddl_p2, nddl_2nd
        integer(kind=8), intent(in) :: nno, nnos
        integer(kind=8), intent(in) :: dimuel, dimdef, kpi
        integer(kind=8), intent(in) :: addeme, addete, addep1, addep2, adde2nd
        real(kind=8), intent(in) :: elem_coor(ndim, nno)
        integer(kind=8), intent(in) :: jv_poids, jv_poids2
        integer(kind=8), intent(in) :: jv_func, jv_func2
        integer(kind=8), intent(in) :: jv_dfunc, jv_dfunc2
        real(kind=8), intent(out) :: dfdi(nno, 3), dfdi2(nnos, 3)
        real(kind=8), intent(out) :: poids, poids2
        real(kind=8), intent(out) :: b(dimdef, dimuel)
    end subroutine cabthm
end interface
