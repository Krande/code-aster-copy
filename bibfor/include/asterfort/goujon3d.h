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
interface 
    subroutine  goujon3d(endo,nbrenf,numr,numf,vecr,&
                           deq,rhor,wpl3,vwpl33,vwpl33t,&
                           rc,sigrf,sigrfissp)
        aster_logical :: endo
        integer(kind=8) :: nbrenf,numr,numf
        real(kind=8) :: vecr(nbrenf,3),deq(nbrenf),rhor(nbrenf)
        real(kind=8) :: wpl3(3),vwpl33(3,3),vwpl33t(3,3)
        real(kind=8), intent(in) :: rc
        real(kind=8) :: sigrf
        real(kind=8) :: sigrfissp(nbrenf,3,3)
    end subroutine goujon3d
end interface 
