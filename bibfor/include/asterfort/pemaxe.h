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
    subroutine pemaxe(resu, nomcha, lieu, nomlie, modele,&
                      chpost, nbcmp, nomcmp, nuord, inst,&
                    nbmail, numemail)
        integer(kind=8) :: nbcmp
        character(len=19) :: resu
        character(len=24) :: nomcha
        character(len=8) :: lieu
        character(len=24) :: nomlie
        character(len=8) :: modele
        character(len=19) :: chpost
        character(len=8) :: nomcmp(nbcmp)
        integer(kind=8) :: nuord
        real(kind=8) :: inst
        integer(kind=8) :: nbmail, numemail(*)
    end subroutine pemaxe
end interface
