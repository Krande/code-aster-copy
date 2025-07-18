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
#include "asterf_types.h"
!
interface
subroutine mdallo(nomres, typcal, nbsauv, base, nbmodes,&
                  rigi, mass, amor, jordr, jdisc,&
                  nbsym, nomsym, jdepl, jvite, jacce,&
                  method, dt, jptem, nbnli, sd_nl_,&
                  jvint, sauve, sd_index, checkarg)
        character(len=8) , intent(in) :: nomres
        character(len=4) , intent(in) :: typcal
        integer(kind=8)          , intent(in) :: nbsauv
        character(len=*) , optional, intent(in)  :: base
        integer(kind=8)          , optional, intent(in)  :: nbmodes
        character(len=*) , optional, intent(in)  :: rigi, mass, amor
        integer(kind=8)          , optional, intent(out) :: jordr, jdisc
        integer(kind=8)          , optional, intent(in)  :: nbsym
        character(len=4) , optional, intent(in)  :: nomsym(*)
        integer(kind=8)          , optional, intent(out) :: jdepl, jvite, jacce
        character(len=*) , optional, intent(in)  :: method
        real(kind=8)     , optional, intent(in)  :: dt
        integer(kind=8)          , optional, intent(out) :: jptem
        integer(kind=8)          , optional, intent(in)  :: nbnli
        character(len=*) , optional, intent(in)  :: sd_nl_
        integer(kind=8)          , optional, intent(out) :: jvint
        character(len=4) , optional, intent(in)  :: sauve
        integer(kind=8)          , optional, intent(in)  :: sd_index
        aster_logical    , optional, intent(in)  :: checkarg
    end subroutine mdallo
end interface
