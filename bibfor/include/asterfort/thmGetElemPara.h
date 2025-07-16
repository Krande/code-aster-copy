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
    subroutine thmGetElemPara(ds_thm   , l_axi    , &
                              type_elem, inte_type, ndim     ,&
                              mecani   , press1   , press2   , tempe  , second, &
                              dimdep   , dimdef   , dimcon   , dimuel ,&
                              nddls    , nddlm    , &
                              nddl_meca, nddl_p1, nddl_p2, nddl_2nd,&
                              nno      , nnos     , &
                              npi      , npg      ,&
                              jv_poids , jv_func  , jv_dfunc ,&
                              jv_poids2, jv_func2 , jv_dfunc2,&
                              jv_gano)
        use THM_type
        type(THM_DS), intent(inout) :: ds_thm
        aster_logical, intent(out) :: l_axi
        character(len=8), intent(out) :: type_elem(2)
        character(len=3), intent(out) :: inte_type
        integer(kind=8), intent(out) :: ndim
        integer(kind=8), intent(out) :: mecani(5), press1(7), press2(7), tempe(5), second(5)
        integer(kind=8), intent(out) :: dimdep, dimdef, dimcon, dimuel
        integer(kind=8), intent(out) :: nddls, nddlm, nddl_meca, nddl_p1, nddl_p2, nddl_2nd
        integer(kind=8), intent(out) :: nno, nnos
        integer(kind=8), intent(out) :: npi, npg
        integer(kind=8), intent(out) :: jv_func, jv_dfunc, jv_poids
        integer(kind=8), intent(out) :: jv_func2, jv_dfunc2, jv_poids2
        integer(kind=8), intent(out) :: jv_gano
    end subroutine thmGetElemPara
end interface
