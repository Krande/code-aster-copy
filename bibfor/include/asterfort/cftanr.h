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
    subroutine cftanr(noma, ndimg, ds_contact, izone,&
                      posnoe, typenm, posenm, numenm, ksipr1,&
                      ksipr2, tau1m, tau2m, tau1, tau2)
        use NonLin_Datastructure_type
        character(len=8) :: noma
        integer(kind=8) :: ndimg
        type(NL_DS_Contact), intent(in) :: ds_contact
        integer(kind=8) :: izone
        integer(kind=8) :: posnoe
        character(len=4) :: typenm
        integer(kind=8) :: posenm
        integer(kind=8) :: numenm
        real(kind=8) :: ksipr1
        real(kind=8) :: ksipr2
        real(kind=8) :: tau1m(3)
        real(kind=8) :: tau2m(3)
        real(kind=8) :: tau1(3)
        real(kind=8) :: tau2(3)
    end subroutine cftanr
end interface
