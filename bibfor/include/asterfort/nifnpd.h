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
    subroutine nifnpd(ndim, nno1, nno2, nno3, npg,&
                      iw, vff1, vff2, vff3, idff1,&
                      vu, vg, vp, typmod, geomi,&
                      sig, ddl, vect)
        integer(kind=8) :: npg
        integer(kind=8) :: nno3
        integer(kind=8) :: nno2
        integer(kind=8) :: nno1
        integer(kind=8) :: ndim
        integer(kind=8) :: iw
        real(kind=8) :: vff1(nno1, npg)
        real(kind=8) :: vff2(nno2, npg)
        real(kind=8) :: vff3(nno3, npg)
        integer(kind=8) :: idff1
        integer(kind=8) :: vu(3, 27)
        integer(kind=8) :: vg(27)
        integer(kind=8) :: vp(27)
        character(len=8) :: typmod(*)
        real(kind=8) :: geomi(ndim, nno1)
        real(kind=8) :: sig(2*ndim+1, npg)
        real(kind=8) :: ddl(*)
        real(kind=8) :: vect(*)
    end subroutine nifnpd
end interface
