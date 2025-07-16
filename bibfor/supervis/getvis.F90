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

subroutine getvis(motfac, motcle, iocc, nbval, vect, &
                  scal, nbret)
! person_in_charge: mathieu.courtois at edf.fr
    implicit none
    character(len=*), intent(in) :: motfac
    character(len=*), intent(in) :: motcle
    integer(kind=8), intent(in), optional :: iocc
    integer(kind=8), intent(in), optional :: nbval
    integer(kind=8), intent(inout), optional :: vect(*)
    integer(kind=8), intent(inout), optional :: scal
    integer(kind=8), intent(out), optional :: nbret
#include "asterc/getvis_wrap.h"
#include "asterfort/assert.h"
!
!   really used variables
    integer(kind=8) :: uioc, unbret, umax
    integer(kind=8) :: uvect(1)
    integer(kind=8) :: vdummy(1)
!
!   motfac + iocc
    if (present(iocc)) then
        uioc = iocc
    else
        uioc = 0
    end if
    ASSERT(motfac == ' ' .or. uioc > 0)
!   vect + nbval
    ASSERT(AU_MOINS_UN3(nbret, scal, vect))
    ASSERT(EXCLUS2(vect, scal))
    if (present(nbval)) then
        umax = nbval
    else
        umax = 1
    end if
!
    if (present(vect)) then
        call getvis_wrap(motfac, motcle, uioc, umax, vect, unbret)
    else if (present(scal)) then
        uvect(1) = scal
        call getvis_wrap(motfac, motcle, uioc, umax, uvect, unbret)
        if (unbret .ne. 0) then
            scal = uvect(1)
        end if
    else
        call getvis_wrap(motfac, motcle, uioc, umax, vdummy, unbret)
    end if
!   if the ".capy" can not ensure that at least 'umax' are provided, you must check
!   the number of values really read using the 'nbret' argument
    ASSERT(present(nbret) .or. umax .eq. unbret)
!
    if (present(nbret)) then
        nbret = unbret
    end if
!
end subroutine getvis
