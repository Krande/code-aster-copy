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
    subroutine xrmes2(ndim, nbnase, cpt, in, ivois,&
                      jsigse, nno, nbcmp, jcnset, dsg11,&
                      dsg22, dsg12)
        integer(kind=8) :: nbnase
        integer(kind=8) :: ndim
        integer(kind=8) :: cpt
        integer(kind=8) :: in
        integer(kind=8) :: ivois
        integer(kind=8) :: jsigse
        integer(kind=8) :: nno
        integer(kind=8) :: nbcmp
        integer(kind=8) :: jcnset
        real(kind=8) :: dsg11(nbnase)
        real(kind=8) :: dsg22(nbnase)
        real(kind=8) :: dsg12(nbnase)
    end subroutine xrmes2
end interface
