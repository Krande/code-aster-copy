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
    subroutine cfparz(ds_contact, iliai, coefff, coefpn, coefpt,&
                      coefte, dissup, izone, ip, numnoe,&
                      posnoe)
        use NonLin_Datastructure_type
        type(NL_DS_Contact), intent(in) :: ds_contact
        integer(kind=8) :: iliai
        real(kind=8) :: coefff
        real(kind=8) :: coefpn
        real(kind=8) :: coefpt
        real(kind=8) :: coefte
        real(kind=8) :: dissup
        integer(kind=8) :: izone
        integer(kind=8) :: ip
        integer(kind=8) :: numnoe
        integer(kind=8) :: posnoe
    end subroutine cfparz
end interface
