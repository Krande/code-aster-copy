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
    subroutine rectfc(nbmode, nbvect, omeshi, npivot, nblagr,&
                      valpro, nvpro, resufi, resufr, nfreq)
        integer(kind=8) :: nfreq
        integer(kind=8) :: nvpro
        integer(kind=8) :: nbmode
        integer(kind=8) :: nbvect
        complex(kind=8) :: omeshi
        integer(kind=8) :: npivot
        integer(kind=8) :: nblagr
        complex(kind=8) :: valpro(nvpro)
        integer(kind=8) :: resufi(nfreq, *)
        real(kind=8) :: resufr(nfreq, *)
    end subroutine rectfc
end interface
