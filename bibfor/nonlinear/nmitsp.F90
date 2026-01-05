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
!
subroutine nmitsp(ds_print, sddisc, iterNewt, retsup)
!
    use NonLin_Datastructure_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/affich.h"
#include "asterfort/assert.h"
#include "asterfort/nmacex.h"
#include "asterfort/nmimpx.h"
#include "asterfort/nmecrr.h"
#include "asterfort/nmlerr.h"
#include "asterfort/utmess.h"
!
    type(NL_DS_Print), intent(in) :: ds_print
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: iterNewt
    integer(kind=8), intent(out) :: retsup
!
! --------------------------------------------------------------------------------------------------
!
! Non-linear algorithm - Discretization management
!
! Action: add more iterations
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_print        : datastructure for printing parameters
! In  sddisc          : datastructure for time discretization
! In  iterNewt        : index of current Newton iteration
! OUT RETSUP : CODE RETOUR
!               0 ON NE PEUT AJOUTER D'ITERATIONS
!               1 ON FAIT DONC DES ITERATIONS EN PLUS
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lExtrapol
    real(kind=8) :: extrapolVale(4), nbIterSuppReal
    integer(kind=8) :: iterSupp, nbIterSupp, vali(2), nbIter, minIter
!
! --------------------------------------------------------------------------------------------------
!
    retsup = 1
    nbIterSupp = 0
    nbIterSuppReal = 0.d0
    lExtrapol = .false.

! - LECTURE DES INFOS SUR LES CONVERGENCES
    call nmlerr(sddisc, 'MNITER', paraValeI_=minIter)
    call nmlerr(sddisc, 'NBITER', paraValeI_=nbIter)

! - AFFICHAGE
    if (iterNewt .ge. nbIter) then
        retsup = 0
        call utmess('I', 'ITERSUPP_2', si=nbIter)
        goto 999
    end if

! - EXTRAPOLATION LINEAIRE DES RESIDUS
    call nmacex(sddisc, iterNewt, lExtrapol, extrapolVale)
    call nmacex(sddisc, iterNewt, lExtrapol, extrapolVale)

! - CALCUL DE LA CIBLE (NOMBRE D'ITERATIONS)
    if (lExtrapol) then
        nbIterSuppReal = (extrapolVale(1)+extrapolVale(2)*log(extrapolVale(4)))/extrapolVale(3)
        call utmess('I', 'EXTRAPOLATION_11')
    else
        nbIterSuppReal = 0.d0
        retsup = 0
        goto 999
    end if

! - NOMBRE D'ITERATIONS SUPPLEMENTAIRES
    nbIterSupp = nint(nbIterSuppReal)
    call utmess('I', 'ITERSUPP_3', si=nbIterSupp)

! - L'EXTRAPOLATION DONNE UN NOMBRE D'ITERATION < LIMITE ITERATION
    if ((nbIterSuppReal*1.20d0) .lt. minIter) then
        retsup = 0
        call utmess('I', 'ITERSUPP_4')
        goto 999
    end if

! - L'EXTRAPOLATION DONNE UN NOMBRE D'ITERATION > LIMITE ITERATION
    if (nbIterSupp .ge. nbIter) then
        retsup = 0
        vali(1) = nbIterSupp
        vali(2) = nbIter
        call utmess('I', 'ITERSUPP_5', ni=2, vali=vali)
    end if
!
999 continue
!
    if (retsup .eq. 1) then
        call utmess('I', 'ITERSUPP_7')
        call affich('MESSAGE', ' ')
        call nmimpx(ds_print)
        iterSupp = 1
    else if (retsup .eq. 0) then
        call utmess('I', 'ITERSUPP_6')
        iterSupp = 0
    else
        ASSERT(ASTER_FALSE)
    end if
    call nmecrr(sddisc, 'ITERSUP', paraValeI_=iterSupp)
!
end subroutine
