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
#include "asterf_types.h"
!
interface
    subroutine fluimp(itypfl, nivpar, nivdef, melflu, typflu,&
                      nuor, freq, freqi, nbm, vite,&
                      npv, carac, calcul, amoc)
        integer(kind=8) :: npv
        integer(kind=8) :: nbm
        integer(kind=8) :: itypfl
        integer(kind=8) :: nivpar
        integer(kind=8) :: nivdef
        character(len=19) :: melflu
        character(len=8) :: typflu
        integer(kind=8) :: nuor(nbm)
        real(kind=8) :: freq(2*nbm*npv)
        real(kind=8) :: freqi(*)
        real(kind=8) :: vite(npv)
        real(kind=8) :: carac(2)
        aster_logical :: calcul(2)
        real(kind=8) :: amoc(*)
    end subroutine fluimp
end interface
