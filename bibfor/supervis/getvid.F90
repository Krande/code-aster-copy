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

subroutine getvid(motfac, motcle, iocc, nbval, vect, &
                  scal, nbret)
! person_in_charge: mathieu.courtois at edf.fr
    implicit none
    character(len=*), intent(in) :: motfac
    character(len=*), intent(in) :: motcle
    integer(kind=8), intent(in), optional :: iocc
    integer(kind=8), intent(in), optional :: nbval
    character(len=*), intent(inout), optional :: vect(*)
    character(len=*), intent(inout), optional :: scal
    integer(kind=8), intent(out), optional :: nbret
#include "asterc/getvid_wrap.h"
#include "asterfort/assert.h"
!
!   really used variables
    integer(kind=8) :: uioc, unbret, umax
!   len=8 should be sufficient, but the same as in getvtx
    integer(kind=8), parameter :: maxlen = 255
    character(len=maxlen) :: uvect(1)
    character(len=1) :: vdummy(1)
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
        call getvid_wrap(motfac, motcle, uioc, umax, vect, unbret)
    else if (present(scal)) then
        ASSERT(len(scal) .le. maxlen)
        uvect(1) = scal
        call getvid_wrap(motfac, motcle, uioc, umax, uvect, unbret)
        if (unbret .ne. 0) then
            scal = uvect(1) (1:len(scal))
        end if
    else
        call getvid_wrap(motfac, motcle, uioc, umax, vdummy, unbret)
    end if
!   if the ".capy" can not ensure that at least 'umax' are provided, you must check
!   the number of values really read using the 'nbret' argument
    ASSERT(present(nbret) .or. umax .eq. unbret)
!
    if (present(nbret)) then
        nbret = unbret
    end if
!
end subroutine getvid
