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
    subroutine mpmod3(basemo, nommes, nbmesu, nbmtot, vcham,&
                      vnoeud, vrange, vorien, nnoema, ncmpma)
        character(len=8) :: basemo
        character(len=8) :: nommes
        integer(kind=8) :: nbmesu
        integer(kind=8) :: nbmtot
        character(len=24) :: vcham
        character(len=24) :: vnoeud
        character(len=24) :: vrange
        character(len=24) :: vorien
        integer(kind=8) :: nnoema
        integer(kind=8) :: ncmpma
    end subroutine mpmod3
end interface
