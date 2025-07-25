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
    subroutine xstudo(ndime, ninter, npts, nptm, ainter,&
                  nbpi, ip1, ip2, pm1a,pm1b, pm2)
        integer(kind=8) :: ndime  
        integer(kind=8) :: ninter
        integer(kind=8) :: npts
        integer(kind=8) :: nptm
        integer(kind=8) :: nbpi
        integer(kind=8) :: ip1(4)
        integer(kind=8) :: ip2(4)
        integer(kind=8) :: pm1a(4)
        integer(kind=8) :: pm1b(4)
        integer(kind=8) :: pm2(4)
        real(kind=8) :: ainter(*)
    end subroutine xstudo
end interface 
