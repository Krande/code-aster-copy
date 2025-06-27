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
#include "asterf_types.h"
!
interface
    subroutine irgmcg(chamsy, partie, ifi, nomcon, ordr,&
                      nbordr, coord, connx, point, nobj,&
                      nbel, nbcmpi, nomcmp, lresu, para,&
                      nomaou, versio)
        character(len=*) :: chamsy
        character(len=*) :: partie
        integer(kind=8) :: ifi
        character(len=*) :: nomcon
        integer(kind=8) :: ordr(*)
        integer(kind=8) :: nbordr
        real(kind=8) :: coord(*)
        integer(kind=8) :: connx(*)
        integer(kind=8) :: point(*)
        character(len=24) :: nobj(28)
        integer(kind=8) :: nbel(28)
        integer(kind=8) :: nbcmpi
        character(len=*) :: nomcmp(*)
        aster_logical :: lresu
        real(kind=8) :: para(*)
        character(len=8) :: nomaou
        integer(kind=8) :: versio
    end subroutine irgmcg
end interface
