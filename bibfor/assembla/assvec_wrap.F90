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
subroutine assvec_wrap(jvBase, vectAsseZ, &
                       nbVectElem, listVectElem, coefVectElem, &
                       numeDofZ, vectScalType)
!
    implicit none
!
#include "asterfort/assvec.h"
!
    character(len=1), intent(in) :: jvBase
    character(len=*), intent(in) :: vectAsseZ
    integer(kind=8), intent(in) :: nbVectElem
    character(len=*), intent(in) :: listVectElem(nbVectElem)
    real(kind=8), intent(in) :: coefVectElem(nbVectElem)
    character(len=*), intent(in) :: numeDofZ
    integer(kind=8), intent(in) :: vectScalType
!
    call assvec(jvBase, vectAsseZ, &
                nbVectElem, listVectElem, coefVectElem, &
                numeDofZ, &
                vectScalType_=vectScalType)
!
end subroutine
