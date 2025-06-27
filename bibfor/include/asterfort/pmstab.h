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
    subroutine pmstab(sigm, sigp, epsm, deps, nbvari,&
                      vim, vip, iforta, instam, instap,&
                      iter, nbpar, nompar, table, vr,&
                      igrad, valimp, imptgt, dsidep, nomvi,&
                      nbvita)
        real(kind=8) :: sigm(6)
        real(kind=8) :: sigp(6)
        real(kind=8) :: epsm(9)
        real(kind=8) :: deps(9)
        integer(kind=8) :: nbvari
        real(kind=8) :: vim(*)
        real(kind=8) :: vip(*)
        integer(kind=8) :: iforta
        real(kind=8) :: instam
        real(kind=8) :: instap
        integer(kind=8) :: iter
        integer(kind=8) :: nbpar
        character(len=16) :: nompar(*)
        character(len=8) :: table
        real(kind=8) :: vr(*)
        integer(kind=8) :: igrad
        real(kind=8) :: valimp(9)
        integer(kind=8) :: imptgt
        real(kind=8) :: dsidep(*)
        character(len=8) :: nomvi(*)
        integer(kind=8) :: nbvita
    end subroutine pmstab
end interface
