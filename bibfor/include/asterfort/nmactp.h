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
!
interface
    subroutine nmactp(ds_print, sddisc, sderro, ds_contact,&
                      ds_conv , nbiter, numins)
        use NonLin_Datastructure_type
        type(NL_DS_Print), intent(in) :: ds_print
        character(len=24), intent(in) :: sderro
        type(NL_DS_Contact), intent(in) :: ds_contact
        character(len=19), intent(in) :: sddisc
        type(NL_DS_Conv), intent(in) :: ds_conv
        integer(kind=8), intent(in) :: nbiter
        integer(kind=8), intent(in) :: numins
    end subroutine nmactp
end interface
