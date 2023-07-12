! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine metaSteelTRCGetParameters(jv_mater, metaSteelPara)
!
use Metallurgy_type
!
implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/rcadma.h"
#include "asterfort/utmess.h"
!
integer, intent(in) :: jv_mater
type(META_SteelParameters), intent(out) :: metaSteelPara
!
! --------------------------------------------------------------------------------------------------
!
! METALLURGY - Steel
!
! Get parameters for TRC curves
!
! --------------------------------------------------------------------------------------------------
!
! In  jv_mater            : coded material address
! Out metaSteelPara       : material parameters for metallurgy of steel
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: toleTemp = 10.
    integer :: jv_pftrc
    integer :: nbcb1, nbcb2, nblexp
    integer :: icodre, nb_trc, nbHist, nbExp, iHist, iExp
    integer :: jftrc, jtrc
    integer :: iadexp, iadckm, iadtrc, shift
    real(kind=8) :: tempAR3FromMate, tempAR3FromTRC
    real(kind=8) :: tempPrev, tempCurr
    aster_logical :: lCooling
!
! --------------------------------------------------------------------------------------------------
!
    call jevech('PFTRC', 'L', jv_pftrc)
    jftrc   = zi(jv_pftrc)
    jtrc    = zi(jv_pftrc+1)
    call rcadma(jv_mater, 'META_ACIER', 'TRC', iadtrc, icodre, 1)
    nbcb1   = nint(zr(iadtrc+1))
    nbHist = nint(zr(iadtrc+2))
    nbcb2 = nint(zr(iadtrc+1+2+nbcb1*nbHist))
    nblexp = nint(zr(iadtrc+1+2+nbcb1*nbHist+1))
    nb_trc = nint(zr(iadtrc+1+2+nbcb1*nbHist+2+nbcb2*nblexp+1))
    ASSERT(nb_trc .eq. 1)
    iadexp = 5+nbcb1*nbHist
    iadckm = 7+nbcb1*nbHist+nbcb2*nblexp
    metaSteelPara%trc%jv_ftrc = jftrc
    metaSteelPara%trc%jv_trc  = jtrc
    metaSteelPara%trc%iadexp  = iadexp
    metaSteelPara%trc%iadtrc  = iadtrc
    metaSteelPara%trc%nbHist = nbHist

! - Parameters for martensite evolution
    metaSteelPara%trc%martensiteLaw%austeniteMin = zr(iadtrc+iadckm-1+1)
    metaSteelPara%trc%martensiteLaw%akm = zr(iadtrc+iadckm-1+2)
    metaSteelPara%trc%martensiteLaw%bkm = zr(iadtrc+iadckm-1+3)
    metaSteelPara%trc%martensiteLaw%lowerSpeed = zr(iadtrc+iadckm-1+4)

! - Parameters for size of austenite grain
    metaSteelPara%trc%austeniteGrain%dref = zr(iadtrc+iadckm-1+5)
    metaSteelPara%trc%austeniteGrain%a = zr(iadtrc+iadckm-1+6)

! - Check consistency of temperature
    tempAR3FromMate = metaSteelPara%ar3
    shift = 0
    do iHist = 1, nbHist
        nbExp = nint(zr(iadtrc+11+9*(iHist-1)))
        lCooling = ASTER_FALSE
        do iExp = 1, nbExp-1
            tempPrev = zr(iadtrc+iadexp-1+4*(shift+iExp))
            tempCurr = zr(iadtrc+iadexp-1+4*(shift+iExp+1))
            if (iExp .ge. 2) then
                if (tempCurr .le. tempPrev) then
                    lCooling = ASTER_TRUE
                    exit
                end if
            end if
        end do
        if (lCooling) then
            iExp = 1
            tempAR3FromTRC = zr(iadtrc+iadexp-1+4*(shift+iExp))
            if (abs(tempAR3FromMate-tempAR3FromTRC) .gt. toleTemp) then
                call utmess('A', "META1_50", sr=toleTemp)
            end if
        end if
        shift = shift+nbExp
    end do
!
end subroutine
