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
      subroutine dflufin3d(sige6,bw,pw,bg,pg,&
               dsw6,delta,rc,&
               xflu,dfin,cmp1,dfmx2)
        real(kind=8), intent(in) :: sige6(6)
        real(kind=8), intent(in) :: bw
        real(kind=8), intent(in) :: pw
        real(kind=8), intent(in) :: bg
        real(kind=8), intent(in) :: pg
        real(kind=8), intent(in) :: dsw6(6)
        real(kind=8), intent(in) :: delta
        real(kind=8), intent(in) :: rc
        real(kind=8), intent(in) :: xflu
        real(kind=8), intent(out) :: dfin
        real(kind=8), intent(out) :: cmp1
        real(kind=8), intent(in) :: dfmx2
    end subroutine dflufin3d
end interface 
