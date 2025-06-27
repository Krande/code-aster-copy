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
    subroutine lcmaei(fami, kpg, ksp, poum, nmater,&
                      imat, necris, necoul, nbval, valres,&
                      nmat, itbint, nfs, nsg, hsri,&
                      ifa, nomfam, nbsys)
        integer(kind=8) :: nsg
        integer(kind=8) :: nmat
        character(len=*) :: fami
        integer(kind=8) :: kpg
        integer(kind=8) :: ksp
        character(len=*) :: poum
        character(len=16) :: nmater
        integer(kind=8) :: imat
        character(len=16) :: necris
        character(len=16) :: necoul
        integer(kind=8) :: nbval
        real(kind=8) :: valres(nmat)
        integer(kind=8) :: itbint
        integer(kind=8) :: nfs
        real(kind=8) :: hsri(nsg, nsg)
        integer(kind=8) :: ifa
        character(len=16) :: nomfam
        integer(kind=8) :: nbsys
    end subroutine lcmaei
end interface
