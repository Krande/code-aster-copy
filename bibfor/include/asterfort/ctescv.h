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
    subroutine ctescv(cvec1, cvec0, cvec01, cvec00, ndim,&
                      xer)
        integer(kind=8) :: ndim
        complex(kind=8) :: cvec1(ndim)
        complex(kind=8) :: cvec0(ndim)
        complex(kind=8) :: cvec01(ndim)
        complex(kind=8) :: cvec00(ndim)
        real(kind=8) :: xer
    end subroutine ctescv
end interface
