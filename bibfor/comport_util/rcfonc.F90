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

subroutine rcfonc(quest, ktrac, jprol, jvale, nbvale, &
                  sigy, e, nu, p, rp, &
                  rprim, airerp, sieleq, dp)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
!
!
    character(len=1), intent(in) :: quest
    integer(kind=8), intent(in) :: ktrac
    integer(kind=8), intent(in) :: jprol
    integer(kind=8), intent(in) :: jvale
    integer(kind=8), intent(in) :: nbvale
    real(kind=8), optional, intent(in) :: e
    real(kind=8), optional, intent(in) :: nu
    real(kind=8), optional, intent(in) :: sieleq
    real(kind=8), optional, intent(in) :: p
    real(kind=8), optional, intent(out) :: sigy
    real(kind=8), optional, intent(out) :: rp
    real(kind=8), optional, intent(out) :: rprim
    real(kind=8), optional, intent(out) :: airerp
    real(kind=8), optional, intent(out) :: dp
!
! --------------------------------------------------------------------------------------------------
!
! Comportment - Utility
!
! Get information by interpolation on traction curve R(p)
!
! --------------------------------------------------------------------------------------------------
!
! In  Quest   : type of information to get
!               'S' - Elastic yield
!               'V' - R(p), dR(p), A[R(p),p]
!               'E' - R(p), dR(p)/dp, dp, A[R(p),p]
! In  ktrac   : type of traction curve
!                1  - 'TRACTION'
!                2  - 'META_TRACTION'
! In  jprol   : adress of .PROL object for R(p) function
! In  jvale   : adress of .VALE object for R(p) function
! In  nbvale  : number of values in R(p) function
! Out sigy    : elastic yield
! In  e       : Young modulus
! In  nu      : Poisson coefficient
! In  sieleq  : elastic stress
! In  p       : internal variable (plastic strain)
! Out rp      : value of R(p)
! Out rprim   : value of dR(p)/dp at p
! Out airrep  : area under R(p) curve
! Out dp      : cumulutated plastic strain
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: tessup
    character(len=1) :: type_prol
    character(len=24) :: func_name
    integer(kind=8) :: jp, jr, i, i0
    real(kind=8) :: troimu, p0, rp0, pp, equ, airerp_int, rprim_int, rp_int
!
! --------------------------------------------------------------------------------------------------
!
    jp = jvale
    jr = jvale+nbvale
    type_prol = zk24(jprol+4) (2:2)
    func_name = zk24(jprol+5)
!
! - Elastic yield
!
    if (quest .eq. 'S') then
        ASSERT(ktrac .eq. 1)
        sigy = zr(jr)
        goto 999
    end if
!
! - Check
!
    if (p .lt. 0) then
        call utmess('F', 'COMPOR5_59')
    end if
!
! DISTINCTION ENTRE TRACTION ET META_TRACTION CAR
! POUR TRACTION:C EST LA COURBE SIGMA(P)
! LE PREMIER POINT CORRESPOND A LA LIMITE D ELASTICITE
! POUR META_TRACTION:C EST LA COURBE R(P)
! DONC IL Y A UN COUPLE DE MOINS QUE DANS TRACTION
!
    if (ktrac .eq. 1) then
        airerp_int = 0.d0
        tessup = .false.
        do i = 1, nbvale-1
            if (p .lt. zr(jp+i)) then
                i0 = i-1
                goto 20
            end if
            airerp_int = airerp_int+(zr(jr+i)+zr(jr+i-1))*(zr(jp+i)-zr(jp+i-1))
        end do
        tessup = .true.
        if (type_prol .eq. 'E') then
            call utmess('F', 'COMPOR5_60', sk=func_name, sr=zr(jp-1+nbvale))
        end if
        i0 = nbvale-1
20      continue
!
! - CALCUL DES VALEURS DE R(P), R'(P) ET AIRE(P)
!
        if (quest .eq. 'V') then
            if (tessup) then
                if (type_prol .eq. 'L') then
                    rprim_int = (zr(jr+i0)-zr(jr+i0-1))/(zr(jp+i0)-zr(jp+i0-1))
                else
                    rprim_int = 0.d0
                end if
            else
                rprim_int = (zr(jr+i0+1)-zr(jr+i0))/(zr(jp+i0+1)-zr(jp+i0))
            end if
            p0 = zr(jp+i0)
            rp0 = zr(jr+i0)
            rp_int = rp0+rprim_int*(p-p0)
            airerp_int = (airerp_int+(rp0+rp_int)*(p-p0))*0.5d0
            goto 999
        end if
    end if
    if (ktrac .eq. 2) then
        airerp_int = 0.d0
        tessup = .false.
!
! - PARCOURS JUSQU'A P
        do i = 1, nbvale
            if (p .lt. zr(jp+i-1)) then
                i0 = i-1
                goto 21
            end if
        end do
        tessup = .true.
        if (type_prol .eq. 'E') then
            call utmess('F', 'COMPOR5_60', sk=func_name, sr=zr(jp-1+nbvale))
        end if
        i0 = nbvale-1
21      continue
!
! - CALCUL DES VALEURS DE R(P) ET R'(P)
!
        if (quest .eq. 'V') then
            if (tessup) then
                if (type_prol .eq. 'L') then
                    rprim_int = (zr(jr+i0)-zr(jr+i0-1))/(zr(jp+i0)-zr(jp+i0-1))
                else
                    rprim_int = 0.d0
                end if
                rp_int = zr(jr+i0)+rprim_int*(p-zr(jp+i0))
            else
                if (i0 .eq. 0) then
                    rprim_int = zr(jr)/zr(jp)
                    p0 = 0.d0
                    rp0 = 0.d0
                else
                    rprim_int = (zr(jr+i0)-zr(jr+i0-1))/(zr(jp+i0)-zr(jp+i0-1))
                    p0 = zr(jp+i0-1)
                    rp0 = zr(jr+i0-1)
                end if
                rp_int = rp0+rprim_int*(p-p0)
            end if
            goto 999
        end if
    end if
!
! - RESOLUTION DE L'EQUATION R(P+DP) + 3 MU DP = SIELEQ
!
    troimu = 1.5d0*e/(1+nu)
    do i = i0+1, nbvale-1
        equ = zr(jr+i)+troimu*(zr(jp+i)-p)-sieleq
        if (equ .gt. 0) then
            i0 = i-1
            goto 40
        end if
        airerp_int = airerp_int+(zr(jr+i)+zr(jr+i-1))*(zr(jp+i)-zr(jp+i-1))
    end do
    tessup = .true.
    if (type_prol .eq. 'E') then
        call utmess('F', 'COMPOR5_60', sk=func_name, sr=zr(jp-1+nbvale))
    end if
    i0 = nbvale-1
40  continue
!
! - CALCUL DES VALEURS DE DP, R(P+DP), R'(P+DP) ET AIRE(P+DP)
!
    if (tessup) then
        if (type_prol .eq. 'L') then
            rprim_int = (zr(jr+i0)-zr(jr+i0-1))/(zr(jp+i0)-zr(jp+i0-1))
        else
            rprim_int = 0.d0
        end if
    else
        rprim_int = (zr(jr+i0+1)-zr(jr+i0))/(zr(jp+i0+1)-zr(jp+i0))
    end if
    p0 = zr(jp+i0)
    rp0 = zr(jr+i0)
    dp = (sieleq-rp0-rprim_int*(p-p0))/(troimu+rprim_int)
    pp = p+dp
    rp_int = rp0+rprim_int*(pp-p0)
    airerp_int = (airerp_int+(rp0+rp_int)*(pp-p0))*0.5d0
!
999 continue
    if (present(airerp)) then
        airerp = airerp_int
    end if
    if (present(rprim)) then
        rprim = rprim_int
    end if
    if (present(rp)) then
        rp = rp_int
    end if
end subroutine
