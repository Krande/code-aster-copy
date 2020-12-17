! --------------------------------------------------------------------
! Copyright (C) 1991 - 2020 - EDF R&D - www.code-aster.org
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
    subroutine cgComputeGtheta(cgField, cgTheta, cgStudy, nume_ordre,&
                               depla, chvite, chacce, time, courb, option, puls, &
                               lmoda, gth)
    use calcG_type
        type(CalcG_field), intent(in) :: cgField
        type(CalcG_theta), intent(in) :: cgTheta
        type(CalcG_Study), intent(in) :: cgStudy
        integer           :: nume_ordre
        character(len=24) :: depla
        character(len=24) :: chvite
        character(len=24) :: chacce
        real(kind=8)      :: time
        character(len=24) :: courb  
        character(len=8)  :: option    
        real(kind=8)      :: puls 
        aster_logical     :: lmoda
        real(kind=8)      :: gth(4)
    end subroutine cgComputeGtheta
end interface
