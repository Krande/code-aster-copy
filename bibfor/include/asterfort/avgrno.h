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
    subroutine avgrno(vwork, tdisp, lisnoe, nbnot, nbordr,&
                      nnoini, nbnop, tspaq, nomcri, nomfor,&
                      grdvie, forvie, fordef, nommai, proaxe,&
                      nommap, cnsr, post, resu)
        integer(kind=8) :: nbnop
        integer(kind=8) :: tdisp
        real(kind=8) :: vwork(tdisp)
        integer(kind=8) :: lisnoe(nbnop)
        integer(kind=8) :: nbnot
        integer(kind=8) :: nbordr
        integer(kind=8) :: nnoini
        integer(kind=8) :: tspaq
        character(len=16) :: nomcri
        character(len=16) :: nomfor
        character(len=16) :: grdvie
        character(len=16) :: forvie
        aster_logical :: fordef
        character(len=8) :: nommai
        character(len=16) :: proaxe
        character(len=8) :: nommap
        character(len=19) :: cnsr
        aster_logical :: post
        real(kind=8) :: resu(7)
    end subroutine avgrno
end interface
