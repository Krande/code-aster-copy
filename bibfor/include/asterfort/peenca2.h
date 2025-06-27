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
    subroutine peenca2(champ, long, vr, nbmail, nummai, ligrel, nbgr, ztot,&
                      option, nbproc, rang)
        integer(kind=8) :: long
        character(len=*) :: champ
        real(kind=8) :: vr(long)
        integer(kind=8) :: nbmail
        integer(kind=8) :: nummai(*)
        character(len=19) :: ligrel
        integer(kind=8) :: nbgr
        real(kind=8) :: ztot
        integer(kind=8) :: option
        integer(kind=8) :: nbproc
        integer(kind=8) :: rang
    end subroutine peenca2
end interface
