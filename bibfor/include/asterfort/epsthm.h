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
    subroutine epsthm(ds_thm   , l_axi    , ndim     ,&
                      addeme   , addep1   , addep2  , addete   , adde2nd, &
                      nno      , nnos     ,&
                      dimuel   , dimdef   , nddls   , nddlm    ,&
                      nddl_meca, nddl_p1  , nddl_p2 , nddl_2nd, &
                      npi      , elem_coor, disp    ,&
                      jv_poids , jv_poids2,&
                      jv_func  , jv_func2 , jv_dfunc, jv_dfunc2,&
                      epsm)
        use THM_type
        type(THM_DS), intent(in) :: ds_thm
        aster_logical, intent(in) :: l_axi
        integer(kind=8), intent(in) :: ndim
        integer(kind=8), intent(in) :: addeme, addep1, addep2, addete, adde2nd
        integer(kind=8), intent(in) :: nno, nnos
        integer(kind=8), intent(in) :: dimuel, dimdef
        integer(kind=8), intent(in) :: nddls, nddlm
        integer(kind=8), intent(in) :: nddl_meca, nddl_p1, nddl_p2, nddl_2nd
        integer(kind=8), intent(in) :: npi
        real(kind=8), intent(in) :: elem_coor(ndim, nno)
        real(kind=8), intent(in) :: disp(*)
        integer(kind=8), intent(in) :: jv_poids, jv_poids2
        integer(kind=8), intent(in) :: jv_func, jv_func2, jv_dfunc, jv_dfunc2
        real(kind=8), intent(out) :: epsm(6,27)
    end subroutine epsthm
end interface
