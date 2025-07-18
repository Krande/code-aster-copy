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

subroutine focoli(ipt, coli, interp, x, y, &
                  rvar, resu, ier)
    implicit none
#include "asterc/r8prem.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ipt, ier
    real(kind=8) :: x(*), y(*), rvar, resu
    character(len=1) :: coli
    character(len=24) :: interp
    real(kind=8) :: valr(3)
! ----------------------------------------------------------------------
#define linlin(x0,x1,y1,x2,y2) y1+(x0-x1)*(y2-y1)/(x2-x1)
#define linlog(x0,x1,y1,x2,y2) exp(log(y1)+(x0-x1)*(log(y2)-log(y1)) \
    /(x2-x1))
#define loglog(x0,x1,y1,x2,y2) exp(log(y1)+(log(x0)-log(x1))*(log(y2) \
    -log(y1))/(log(x2)-log(x1)))
#define loglin(x0,x1,y1,x2,y2) y1+(log(x0)-log(x1))*(y2-y1) \
    /(log(x2)-log(x1))
! ----------------------------------------------------------------------
!
!     --- PAS D'INTERPOLATION ---
!
    ier = 0
    if (coli .eq. 'C') then
        resu = y(ipt)
!
!     --- INTERPOLATION ---
!
    else if (coli .eq. 'I') then
        if (interp .eq. 'LIN LIN ') then
            resu = linlin(rvar, x(ipt), y(ipt), x(ipt+1), y(ipt+1))
!
        else if (interp .eq. 'LIN LOG ') then
            if (y(ipt) .lt. r8prem() .or. y(ipt+1) .lt. r8prem()) then
                ier = 250
                valr(1) = rvar
                valr(2) = x(ipt)
                valr(3) = x(ipt+1)
                call utmess('A', 'UTILITAI6_15', nr=3, valr=valr)
            else
                resu = linlog(rvar, x(ipt), y(ipt), x(ipt+1), y(ipt+1))
            end if
!
        else if (interp .eq. 'LOG LOG ') then
            if (x(ipt) .lt. r8prem() .or. x(ipt+1) .lt. r8prem() .or. &
                y(ipt) .lt. r8prem() .or. y(ipt+1) .lt. r8prem()) then
                ier = 250
                valr(1) = rvar
                valr(2) = x(ipt)
                valr(3) = x(ipt+1)
                call utmess('A', 'UTILITAI6_15', nr=3, valr=valr)
            else
                resu = loglog(rvar, x(ipt), y(ipt), x(ipt+1), y(ipt+1))
            end if
!
        else if (interp .eq. 'LOG LIN ') then
            if (x(ipt) .lt. r8prem() .or. x(ipt+1) .lt. r8prem()) then
                ier = 250
                valr(1) = rvar
                valr(2) = x(ipt)
                valr(3) = x(ipt+1)
                call utmess('A', 'UTILITAI6_15', nr=3, valr=valr)
            else
                resu = loglin(rvar, x(ipt), y(ipt), x(ipt+1), y(ipt+1))
            end if
!
        else
            ier = 230
            call utmess('A', 'UTILITAI_84', sk=interp)
        end if
!
!     --- EXTRAPOLATION ---
!
    else if (coli .eq. 'E') then
        resu = linlin(rvar, x(ipt), y(ipt), x(ipt+1), y(ipt+1))
!
    else
        ier = 240
        call utmess('A', 'UTILITAI_85', sk=coli)
    end if
!
end subroutine
