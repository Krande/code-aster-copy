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
    subroutine xstjon(elrefp, ndim, joncno, jlsn, igeom, nfiss, nfisc, fisco, nnops,&
                      txlsn, n, c)
        character(len=8) :: elrefp
        integer(kind=8) :: ndim
        integer(kind=8) :: joncno
        integer(kind=8) :: jlsn
        integer(kind=8) :: igeom
        integer(kind=8) :: nfiss
        integer(kind=8) :: nfisc
        integer(kind=8) :: fisco(*)
        integer(kind=8) :: nnops
        real(kind=8) :: txlsn(28)
        integer(kind=8), intent(in), optional :: n
        real(kind=8), intent(in), optional :: c(ndim)
    end subroutine xstjon
end interface
