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

subroutine w155m3(numa, jce2d, jce2l, jce2v, isp, &
                  nucou, nusec, nufib, posic, posis)
! person_in_charge: jacques.pellet at edf.fr
! ======================================================================
    implicit none
#include "jeveux.h"
!
#include "asterfort/assert.h"
#include "asterfort/cesexi.h"
    integer(kind=8) :: numa, nucou, nusec, nufib, posic, posis, isp
    integer(kind=8) :: jce2l, jce2d, jce2v
    integer(kind=8) :: nbcou, nbsec, nbfib, isec, icou
    integer(kind=8) :: iad1, iad2, iad3, iad4
    character(len=8) :: typma
!
! ----------------------------------------------------------------------
! BUT : DETERMINER NUCOU, NUSEC, NUFIB, POSIC ET POSIS
!       A PARTIR DU NUMERO DE SOUS-POINT ISP
! ----------------------------------------------------------------------
!
!   ON DETERMINE SI LA MAILLE EST DE TYPE COQUE, PMF OU TUYAU :
!
!   CMP1 = COQ_NCOU
    call cesexi('C', jce2d, jce2l, numa, 1, 1, 1, iad1)
!   CMP2 = TUY_NCOU
    call cesexi('C', jce2d, jce2l, numa, 1, 1, 2, iad2)
!   CMP3 = TUY_NSEC
    call cesexi('C', jce2d, jce2l, numa, 1, 1, 3, iad3)
!   CMP4 = NBFIBR
    call cesexi('C', jce2d, jce2l, numa, 1, 1, 4, iad4)
!
    if (iad4 .gt. 0) then
        typma = 'PMF'
        nbfib = zi(jce2v-1+iad4)
    else if (iad2 .gt. 0) then
        typma = 'TUY'
        ASSERT(iad3 .gt. 0)
        nbcou = zi(jce2v-1+iad2)
        nbsec = zi(jce2v-1+iad3)
    else if (iad1 .gt. 0) then
        typma = 'COQ'
        nbcou = zi(jce2v-1+iad1)
    else
        ASSERT(.false.)
    end if
!
    nucou = -999
    nusec = -999
    nufib = -999
    posic = -999
    posis = -999
!
    if (typma .eq. 'PMF') then
!       -------------------------
        ASSERT(isp .le. nbfib)
        nufib = isp
!
    else if (typma .eq. 'COQ') then
!       -------------------------
        ASSERT(isp .le. nbcou*3)
        nucou = (isp+2)/3
        posic = mod(isp+2, 3)-1
!
    else if (typma .eq. 'TUY') then
!       -------------------------
        ASSERT(isp .le. (2*nbcou+1)*(2*nbsec+1))
        icou = (isp-1)/(2*nbsec+1)+1
        isec = isp-(icou-1)*(2*nbsec+1)
        if (icou .eq. 1) then
            nucou = 1
            posic = -1
        else
            nucou = icou/2
            posic = icou-2*nucou
        end if
        if (isec .eq. 1) then
            nusec = 1
            posis = -1
        else
            nusec = isec/2
            posis = isec-2*nusec
        end if
!
    else
        ASSERT(.false.)
    end if
!
end subroutine
