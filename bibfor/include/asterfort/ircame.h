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
    subroutine ircame(ifi, nochmd, chanom, typech, modele,&
                      nbcmp, nomcmp, etiqcp, partie, numpt,&
                      instan, numord, adsk, adsd, adsc,&
                      adsv, adsl, nbenec, lienec, sdcarm,&
                      carael, field_type, nbCmpDyna, lfichUniq, codret)
        integer(kind=8) :: ifi
        character(len=64) :: nochmd
        character(len=19) :: chanom
        character(len=8) :: typech
        character(len=8) :: modele
        integer(kind=8) :: nbcmp
        character(len=*) :: nomcmp(*)
        character(len=*) :: etiqcp
        character(len=*) :: partie
        integer(kind=8) :: numpt
        real(kind=8) :: instan
        integer(kind=8) :: numord
        integer(kind=8) :: adsk
        integer(kind=8) :: adsd
        integer(kind=8) :: adsc
        integer(kind=8) :: adsv
        integer(kind=8) :: adsl
        integer(kind=8) :: nbenec
        integer(kind=8) :: lienec(*)
        character(len=8) :: sdcarm, carael
        character(len=16), intent(in) :: field_type
        integer(kind=8), intent(inout) :: nbCmpDyna
        aster_logical :: lfichUniq
        integer(kind=8) :: codret
    end subroutine ircame
end interface
