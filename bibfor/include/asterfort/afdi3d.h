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
    subroutine afdi3d(irep, eta, car, val, jdc,&
                      jdv, ivr, iv, kma, ncmp,&
                      ntp, jdcinf, jdvinf, isym )
        integer(kind=8) :: irep
        real(kind=8) :: eta
        character(len=*) :: car
        real(kind=8) :: val(*)
        integer(kind=8) :: jdc(3)
        integer(kind=8) :: jdv(3)
        integer(kind=8) :: ivr(*)
        integer(kind=8) :: iv
        character(len=1) :: kma(3)
        integer(kind=8) :: ncmp
        integer(kind=8) :: ntp
        integer(kind=8) :: jdcinf
        integer(kind=8) :: jdvinf
        integer(kind=8) :: isym
    end subroutine afdi3d
end interface
