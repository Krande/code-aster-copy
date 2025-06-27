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
    subroutine compute_ineq_conditions_vector(jsecmb, nbliai, neq   ,&
                                            japptr, japddl, japcoe,&
                                            jjeux , jtacf , jmu   ,&
                                            njeux , ztacf , iliac )
        integer(kind=8) :: jsecmb
        real(kind=8) :: jeuini, coefpn, lambdc
        integer(kind=8) :: iliai, iliac
        integer(kind=8) :: jmu
        integer(kind=8) :: japcoe, japddl, japptr
        integer(kind=8) :: jtacf
        integer(kind=8) :: jjeux, njeux
        integer(kind=8) :: ztacf
        integer(kind=8) :: nbliai, neq
        integer(kind=8) :: nbddl, jdecal
    end subroutine compute_ineq_conditions_vector
end interface
