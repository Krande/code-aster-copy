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
subroutine rcfode(ifon, temp, f, df)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ifon
    real(kind=8) :: temp, f, df
! ......................................................................
!     OBTENTION DE LA VALEUR DE LA FONCTION ET DE SA DERIVEE POUR UNE
!     FONCTION DE LA TEMPERATURE LINEAIRE PAR MORCEAU
! IN   IFON   : ADRESSE DANS LE MATERIAU CODE DE LA FONCTION
! IN   TEMP   : TEMPERATURE AU POINT DE GAUSS CONSIDERE
! OUT  F      : VALEUR DE LA FONCTION
! OUT  DF     : VALEUR DE LA DERIVEE DE LA FONCTION
!
!
    real(kind=8) :: valr
!
!
    integer(kind=8) :: jpro, jvalf, jv, jp, nbvf
    aster_logical :: tesinf, tessup, entre, deja, avant
! ----------------------------------------------------------------------
! PARAMETER ASSOCIE AU MATERIAU CODE
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ideb, ifin, incr, indfct, isave
!-----------------------------------------------------------------------
    parameter(indfct=7)
! DEB ------------------------------------------------------------------
    jpro = zi(ifon+1)
    jvalf = zi(ifon+2)
    if (zk24(jpro) (1:1) .eq. 'C') then
        f = zr(jvalf+1)
        df = 0.d0
        goto 101
    else if (zk24(jpro) (1:1) .eq. 'I') then
        call utmess('F', 'MODELISA6_31')
    else if (zk24(jpro) (1:1) .eq. 'N') then
        call utmess('F', 'MODELISA6_58')
    end if
    nbvf = zi(ifon)
    isave = zi(ifon+indfct)
!
    deja = temp .ge. zr(jvalf+isave-1) .and. temp .le. zr(jvalf+isave)
    avant = temp .lt. zr(jvalf+isave-1)
    tesinf = temp .lt. zr(jvalf)
    tessup = temp .gt. zr(jvalf+nbvf-1)
    entre = .not. tesinf .and. .not. tessup
    if (deja) then
        jp = jvalf+isave
        jv = jp+nbvf
        df = (zr(jv)-zr(jv-1))/(zr(jp)-zr(jp-1))
        f = df*(temp-zr(jp-1))+zr(jv-1)
        goto 100
    else
        if (avant) then
            ifin = jvalf
            ideb = jvalf+isave-1
            incr = -1
        else
            ideb = jvalf+isave
            ifin = jvalf+nbvf-1
            incr = 1
        end if
    end if
    if (entre) then
        if (incr .gt. 0) then
            do jp = ideb, ifin, incr
                jv = jp+nbvf
                if (zr(jp) .ge. temp) then
                    df = (zr(jv)-zr(jv-1))/(zr(jp)-zr(jp-1))
                    f = df*(temp-zr(jp-1))+zr(jv-1)
                    isave = jp-jvalf
                    goto 5
                end if
            end do
5           continue
        else
            do jp = ideb, ifin, incr
                jv = jp+nbvf
                if (zr(jp) .le. temp) then
                    df = (zr(jv+1)-zr(jv))/(zr(jp+1)-zr(jp))
                    f = df*(temp-zr(jp))+zr(jv)
                    isave = jp-jvalf+1
                    goto 6
                end if
            end do
6           continue
        end if
    else if (tesinf) then
        jv = jvalf+nbvf
        jp = jvalf
        if (zk24(jpro+4) (2:2) .eq. 'C') then
            df = 0.0d0
            f = zr(jv)
        else if (zk24(jpro+4) (1:1) .eq. 'L') then
            df = (zr(jv+1)-zr(jv))/(zr(jp+1)-zr(jp))
            f = df*(temp-zr(jp))+zr(jv)
        else if (zk24(jpro+4) (1:1) .eq. 'E') then
            valr = temp
            call utmess('F', 'MODELISA8_93', sr=valr)
        end if
        isave = 1
    else if (tessup) then
        jv = jvalf+2*nbvf-1
        jp = jvalf+nbvf-1
        if (zk24(jpro+4) (2:2) .eq. 'C') then
            df = 0.0d0
            f = zr(jv)
        else if (zk24(jpro+4) (2:2) .eq. 'L') then
            df = (zr(jv)-zr(jv-1))/(zr(jp)-zr(jp-1))
            f = df*(temp-zr(jp-1))+zr(jv-1)
        else if (zk24(jpro+4) (2:2) .eq. 'E') then
            valr = temp
            call utmess('F', 'MODELISA8_94', sr=valr)
        end if
        isave = nbvf-1
    end if
100 continue
    zi(ifon+indfct) = isave
101 continue
! FIN ------------------------------------------------------------------
end subroutine
