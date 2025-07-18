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
subroutine rcfodi(ifon, beta, f, df)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ifon
    real(kind=8) :: beta, f, df
! ......................................................................
!     ROUTINE INVERSE DE RCFODE :
!     UTILISEE UNIQUEMENT POUR LA FONCTION ENTHALPIE DANS
!     L OPERATEUR THER_NON_LINE_MO (REPERE MOBILE)
!     LA FONCTION FOURNIE PAR LE MATERIAU CODE EST BETA=F(TEMPERATURE)
!     ON CALCULE ICI TEMPERATURE=F-1(BETA), DE MEME QUE LA DERIVEE DE
!     CETTE FONCTION.
! IN   IFON   : ADRESSE DANS LE MATERIAU CODE DE LA FONCTION
! IN   BETA   : ENTHALPIE AU POINT DE GAUSS CONSIDERE
! OUT  F      : VALEUR DE LA FONCTION (TEMPERATURE)
! OUT  DF     : VALEUR DE LA DERIVEE DE LA FONCTION
!
!
!
!
    integer(kind=8) :: jpro, jvalf, jv, jp, nbvf
    integer(kind=8) :: isave, ideb, ifin, incr, indfct
    aster_logical :: tesinf, tessup, entre, deja, avant
! ----------------------------------------------------------------------
! PARAMETER ASSOCIE AU MATERIAU CODE
!
    parameter(indfct=7)
! DEB ------------------------------------------------------------------
    nbvf = zi(ifon)
    jpro = zi(ifon+1)
    jvalf = zi(ifon+2)
    if (zk24(jpro) (1:1) .eq. 'C') then
        f = zr(jvalf+nbvf+1)
        df = 0.d0
        goto 101
    else if (zk24(jpro) (1:1) .eq. 'N') then
        call utmess('F', 'MODELISA6_58')
    end if
    isave = zi(ifon+indfct)
!
    deja = beta .ge. zr(jvalf+isave+nbvf-1) .and. beta .le. zr(jvalf+isave+nbvf)
    avant = beta .lt. zr(jvalf+isave+nbvf-1)
    tesinf = beta .lt. zr(jvalf+nbvf)
    tessup = beta .gt. zr(jvalf+nbvf+nbvf-1)
    entre = .not. tesinf .and. .not. tessup
    if (deja) then
        jp = jvalf+isave
        jv = jp+nbvf
        df = (zr(jp)-zr(jp-1))/(zr(jv)-zr(jv-1))
        f = df*(beta-zr(jv-1))+zr(jp-1)
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
                if (zr(jv) .ge. beta) then
                    df = (zr(jp)-zr(jp-1))/(zr(jv)-zr(jv-1))
                    f = df*(beta-zr(jv-1))+zr(jp-1)
                    isave = jp-jvalf
                    goto 5
                end if
            end do
5           continue
        else
            do jp = ideb, ifin, incr
                jv = jp+nbvf
                if (zr(jv) .le. beta) then
                    df = (zr(jp+1)-zr(jp))/(zr(jv+1)-zr(jv))
                    f = df*(beta-zr(jv))+zr(jp)
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
            f = zr(jp)
        else if (zk24(jpro+4) (1:1) .eq. 'L') then
            df = (zr(jp+1)-zr(jp))/(zr(jv+1)-zr(jv))
            f = df*(beta-zr(jv))+zr(jp)
        else if (zk24(jpro+4) (1:1) .eq. 'E') then
            call utmess('F', 'MODELISA4_63')
        end if
        isave = 1
    else if (tessup) then
        jv = jvalf+2*nbvf-1
        jp = jvalf+nbvf-1
        if (zk24(jpro+4) (2:2) .eq. 'C') then
            df = 0.0d0
            f = zr(jp)
        else if (zk24(jpro+4) (2:2) .eq. 'L') then
            df = (zr(jp)-zr(jp-1))/(zr(jv)-zr(jv-1))
            f = df*(beta-zr(jv-1))+zr(jp-1)
        else if (zk24(jpro+4) (2:2) .eq. 'E') then
            call utmess('F', 'MODELISA4_65')
        end if
        isave = nbvf-1
    end if
100 continue
    zi(ifon+indfct) = isave
101 continue
! FIN ------------------------------------------------------------------
end subroutine
