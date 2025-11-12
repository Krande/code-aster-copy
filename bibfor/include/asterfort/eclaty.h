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
interface
    subroutine eclaty(typeElemName, &
                      elrefa, fapg, &
                      mxnbn2, mxnbpi, mxnbte, mxnbse, &
                      npg, npoini, &
                      nterm1, nsomm1, csomm1, &
                      typeCellNume, nbno2, connx, &
                      nbsel, corsel, iret)
        character(len=16), intent(in) :: typeElemName
        character(len=8), intent(in) :: elrefa, fapg
        integer(kind=8), intent(in) :: mxnbn2, mxnbpi, mxnbte, mxnbse
        integer(kind=8), intent(out) :: npg, npoini
        integer(kind=8), intent(out) :: nterm1(mxnbpi), nsomm1(mxnbpi, mxnbte)
        real(kind=8), intent(out) :: csomm1(mxnbpi, mxnbte)
        integer(kind=8), intent(out) :: typeCellNume(mxnbse), nbno2(mxnbse), connx(mxnbn2, mxnbse)
        integer(kind=8), intent(out) :: nbsel, corsel(mxnbse), iret
    end subroutine eclaty
end interface
