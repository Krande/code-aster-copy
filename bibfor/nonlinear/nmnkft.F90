! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
! aslint: disable=W0413
!
subroutine nmnkft(solver, sddisc, iterNewt_)
!
    implicit none
!
#include "asterfort/jeveuo.h"
#include "asterfort/nmlere.h"
#include "asterfort/nmecrr.h"
#include "asterfort/nmlerr.h"
!
    character(len=19), intent(in) :: solver, sddisc
    integer(kind=8), optional, intent(in) :: iterNewt_
!
! --------------------------------------------------------------------------------------------------
!
! Non-linear algorithm - Discretization management
!
! Action: update Newto-Krylov
!
! --------------------------------------------------------------------------------------------------
!
! SCHEMA DE CALCUL INPIRE DE "SOLVING NONLINEAR EQUATION WITH
!     NEWTON'S METHOD", C.T. KELLEY, SIAM, PAGE 62-63
!
! In  solver           : name of datastructure for solver
! In  sddisc           : datastructure for time discretization
! In  iterNewt         : index of current Newton iteration
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iterNewt
    real(kind=8) :: newtKrylResi, newtKrylResiPrev, resiCurr(1), resiPrev(1), solvMini
    real(kind=8), pointer :: slvr(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jeveuo(solver//'.SLVR', 'E', vr=slvr)
    if (.not. present(iterNewt_)) then
        call nmlerr(sddisc, 'INIT_NEWTON_KRYLOV', paraValeR_=newtKrylResi)
    else
        iterNewt = iterNewt_
        if (iterNewt .eq. 0) then
            call nmlere(sddisc, 'L', 'VCHAR', iterNewt, resiPrev(1))
        else
            call nmlere(sddisc, 'L', 'VMAXI', iterNewt-1, resiPrev(1))
        end if
        call nmlerr(sddisc, 'ITER_NEWTON_KRYLOV', paraValeR_=newtKrylResiPrev)
        call nmlere(sddisc, 'L', 'VMAXI', iterNewt, resiCurr(1))
        if (resiPrev(1) .eq. 0.d0) then
            newtKrylResi = newtKrylResiPrev
            goto 10
        end if
        if ((0.9d0*newtKrylResiPrev**2) .gt. 0.2d0) then
            newtKrylResi = &
                min(max(0.1d0*resiCurr(1)**2/resiPrev(1)**2, 0.9d0*newtKrylResiPrev**2), &
                    4.d-1*newtKrylResiPrev)
        else
            solvMini = slvr(1)
            newtKrylResi = &
                max(min(0.1d0*resiCurr(1)**2/resiPrev(1)**2, 4.d-1*newtKrylResiPrev), &
                    solvMini)
        end if
    end if
!
10  continue

! - STOCKAGE DE LA PRECISION CALCULEE POUR ITERATION SUIVANTE
    call nmecrr(sddisc, 'ITER_NEWTON_KRYLOV', paraValeR_=newtKrylResi)

! - COPIE DE LA PRECISION CALCULEE DANS LA SD SOLVEUR
    slvr(2) = newtKrylResi
!
end subroutine
