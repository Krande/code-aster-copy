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
    subroutine nmfext(eta, listFuncActi, veasse, cnfext,&
                      ds_contact_, sddyna_, nlDynaDamping_)
        use NonLin_Datastructure_type
        use NonLinearDyna_type
        real(kind=8), intent(in) :: eta
        integer(kind=8), intent(in) :: listFuncActi(*)
        character(len=19) :: veasse(*)
        type(NL_DS_Contact), optional, intent(in) :: ds_contact_
        character(len=19), intent(in) :: cnfext
        character(len=19), optional, intent(in) :: sddyna_
        type(NLDYNA_DAMPING), optional,intent(in) :: nlDynaDamping_
    end subroutine nmfext
end interface
