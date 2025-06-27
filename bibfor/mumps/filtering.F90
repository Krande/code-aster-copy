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

subroutine filtering(ltypr, k, lmnsy, jvale, jvale2, rfiltre, rmin, rmax, &
                     nfilt1, nfilt2, nfilt3, nzloc, iok)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
!
    aster_logical, intent(in) :: ltypr, lmnsy
    integer(kind=8), intent(in) :: k, jvale, jvale2
    real(kind=8), intent(in) :: rfiltre, rmin, rmax
    integer(kind=8), intent(inout) :: nfilt1, nfilt2, nfilt3, nzloc
    integer(kind=4), intent(inout) :: iok(*)
!
!---------------------------------------------------------------------------------------------------
!
! Le but est de construire le graphe de comm optimisé pour les comm point à point
!
!---------------------------------------------------------------------------------------------------
!
    real(kind=8) :: rtest, rval
    complex(kind=8) :: cval
!
    if (ltypr) then
        rval = zr(jvale-1+k)
        if (lmnsy) rval = zr(jvale2-1+k)
        if (rval /= 0.d0) then
            rtest = abs(rval)
            if (rtest .gt. rfiltre) then
                if (rtest .lt. rmin) then
                    nfilt3 = nfilt3+1
                    iok(k) = -2
                else if (rtest .gt. rmax) then
                    nfilt2 = nfilt2+1
                    iok(k) = -1
                else
                    iok(k) = 1
                end if
                nzloc = nzloc+1
            end if
        else
! ---   TERME RIGOUREUSEMENT NUL
            nfilt1 = nfilt1+1
        end if
!
    else
        cval = zc(jvale-1+k)
        if (lmnsy) cval = zc(jvale2-1+k)
        if (cval /= (0.d0, 0.d0)) then
            rtest = abs(cval)
            if (rtest .gt. rfiltre) then
                if (rtest .lt. rmin) then
                    nfilt3 = nfilt3+1
                    iok(k) = -2
                else if (rtest .gt. rmax) then
                    nfilt2 = nfilt2+1
                    iok(k) = -1
                else
                    iok(k) = 1
                end if
                nzloc = nzloc+1
            end if
        else
! ---   TERME RIGOUREUSEMENT NUL
            nfilt1 = nfilt1+1
        end if
    end if
!
end subroutine
