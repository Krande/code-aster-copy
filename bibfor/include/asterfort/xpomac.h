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
    subroutine xpomac(malini, mailc, listno, nbnoc, nbmac,&
                      maxfem, nivgrm, cns1, cns2, ces1,&
                      ces2, cesvi1, cesvi2, resuco, comps1,&
                      comps2, pre1)
        character(len=8) :: malini
        character(len=24) :: mailc
        character(len=24) :: listno
        integer(kind=8) :: nbnoc
        integer(kind=8) :: nbmac
        character(len=8) :: maxfem
        character(len=24) :: nivgrm
        character(len=19) :: cns1
        character(len=19) :: cns2
        character(len=19) :: ces1
        character(len=19) :: ces2
        character(len=19) :: cesvi1
        character(len=19) :: cesvi2
        character(len=8) :: resuco
        character(len=19) :: comps1
        character(len=19) :: comps2
        aster_logical :: pre1
    end subroutine xpomac
end interface 
