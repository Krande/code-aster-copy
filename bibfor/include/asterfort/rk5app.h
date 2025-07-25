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
#include "asterf_types.h"
!
interface
    subroutine rk5app(nbeq, vparam_real, vparam_int, vparam_car, dtemps, &
                      yinit, dyinit, rk5fct, solu, decoup)
        integer(kind=8)          :: nbeq
        real(kind=8)     :: vparam_real(*)
        integer(kind=8)          :: vparam_int(*)
        character(len=*) :: vparam_car(*)
        real(kind=8)     :: dtemps
        real(kind=8)     :: yinit(nbeq)
        real(kind=8)     :: dyinit(nbeq)
        real(kind=8)     :: solu(3*nbeq)
        aster_logical    :: decoup
!
        interface
            subroutine rk5fct(ppr, ppi, ppc, yy0, dy0, dyy, decoup)
                real(kind=8)     :: ppr(*)
                integer(kind=8)          :: ppi(*)
                character(len=*) :: ppc(*)
                real(kind=8)     :: yy0(*)
                real(kind=8)     :: dy0(*)
                real(kind=8)     :: dyy(*)
                aster_logical    :: decoup
            end subroutine rk5fct
        end interface
!
    end subroutine rk5app
end interface
