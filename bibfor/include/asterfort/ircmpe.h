! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
#include "asterf_types.h"
!
interface
    subroutine ircmpe(nofimd, ncmpve, numcmp, exicmp, nbvato,&
                      nbmaec, limaec, adsd, adsl, nbimpr,&
                      ncaimi, ncaimk, tyefma, typmai, typgeo,&
                      nomtyp, typech, profas, promed, prorec,&
                      nroimp, chanom, sdcarm, field_type, nosdfu)
        integer :: nbvato
        integer :: ncmpve
        character(len=*) :: nofimd
        integer :: numcmp(ncmpve)
        aster_logical :: exicmp(nbvato)
        integer :: nbmaec
        integer :: limaec(*)
        integer :: adsd
        integer :: adsl
        integer :: nbimpr
        character(len=24) :: ncaimi
        character(len=24) :: ncaimk
        integer :: tyefma(*)
        integer :: typmai(*)
        integer :: typgeo(*)
        character(len=8) :: nomtyp(*)
        character(len=8) :: typech
        integer :: profas(nbvato)
        integer :: promed(nbvato)
        integer :: prorec(nbvato)
        integer :: nroimp(nbvato)
        character(len=19) :: chanom
        character(len=8) :: sdcarm
        character(len=16), intent(in) :: field_type
        character(len=8) :: nosdfu
    end subroutine ircmpe
end interface
