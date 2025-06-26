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
#include "asterfort/Behaviour_type.h"
#include "asterf_types.h"
!
interface
    subroutine tgveri(option, carcri, compor, nno, geom, &
                      ndim, nddl, deplp, sdepl, vectu, &
                      svect, ncont, contp, scont, nvari, &
                      varip, svari, matuu, smatr, matsym, &
                      epsilo, varia, iret)
        character(len=16) :: option
        real(kind=8) :: carcri(CARCRI_SIZE)
        character(len=16) :: compor(COMPOR_SIZE)
        integer(kind=8) :: nno
        real(kind=8) :: geom(*)
        integer(kind=8) :: ndim
        integer(kind=8) :: nddl
        real(kind=8) :: deplp(*)
        real(kind=8) :: sdepl(*)
        real(kind=8) :: vectu(*)
        real(kind=8) :: svect(*)
        integer(kind=8) :: ncont
        real(kind=8) :: contp(*)
        real(kind=8) :: scont(*)
        integer(kind=8) :: nvari
        real(kind=8) :: varip(*)
        real(kind=8) :: svari(*)
        real(kind=8) :: matuu(*)
        real(kind=8) :: smatr(*)
        aster_logical :: matsym
        real(kind=8) :: varia(*)
        real(kind=8) :: epsilo
        integer(kind=8) :: iret
    end subroutine tgveri
end interface
