! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
      subroutine conso3d(epsmk,epser,ccmin,ccmax,&
                       epsm6,epse6,cc3,vepsm33,vepsm33t,CWp,CMp,CTD,CTV)
        real(kind=8), intent(in) :: epsmk
        real(kind=8), intent(in) :: epser
        real(kind=8), intent(out) :: ccmin
        real(kind=8), intent(out) :: ccmax
        real(kind=8), intent(in) :: epsm6(6)
        real(kind=8), intent(in) :: epse6(6)
        real(kind=8), intent(out) :: cc3(3)
        real(kind=8), intent(out) :: vepsm33(3,3)
        real(kind=8), intent(out) :: vepsm33t(3,3)
        real(kind=8), intent(in) :: CWp
        real(kind=8), intent(in) :: CMp
        real(kind=8), intent(in) :: CTD
        real(kind=8), intent(in) :: CTV
    end subroutine conso3d
end interface 
