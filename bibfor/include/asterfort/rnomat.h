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
    subroutine rnomat(icesd, icesl, icesv, imap, nomcri,&
                      adrma, jtypma, k, optio, vala,&
                      valb, coefpa, nommat)
        integer(kind=8) :: icesd
        integer(kind=8) :: icesl
        integer(kind=8) :: icesv
        integer(kind=8) :: imap
        character(len=16) :: nomcri
        integer(kind=8) :: adrma
        integer(kind=8) :: jtypma
        integer(kind=8) :: k
        character(len=10) :: optio
        real(kind=8) :: vala
        real(kind=8) :: valb
        real(kind=8) :: coefpa
        character(len=8) :: nommat
    end subroutine rnomat
end interface
