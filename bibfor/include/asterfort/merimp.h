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
interface
    subroutine merimp(l_xfem         , l_dyna     , &
                      model          , cara_elem  , sddyna  , iter_newt,&
                      ds_constitutive, ds_material,&
                      hval_incr      , hval_algo  , caco3d   ,&
                      mxchin         , lpain      , lchin   , nbin)
        use NonLin_Datastructure_type
        aster_logical, intent(in) :: l_xfem, l_dyna
        character(len=24), intent(in) :: model, cara_elem
        character(len=19), intent(in) :: sddyna
        integer(kind=8), intent(in) :: iter_newt
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        type(NL_DS_Material), intent(in) :: ds_material
        character(len=19), intent(in) :: hval_incr(*), hval_algo(*)
        character(len=24), intent(in) :: caco3d
        integer(kind=8), intent(in) :: mxchin
        character(len=8), intent(inout) :: lpain(mxchin)
        character(len=19), intent(inout) :: lchin(mxchin)
        integer(kind=8), intent(out) :: nbin
    end subroutine merimp
end interface
