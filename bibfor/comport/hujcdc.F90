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

subroutine hujcdc(k, mater, sig, vin, seuil)
    implicit none
!    HUJEUX:  SEUIL DU MECANISME DEVIATOIRE CYCLIQUE K(=1 A 3)
!             FD(K) = QIIC(K) + M*PK*RK*( 1 - B*LOG(PK/PC) )
!             QIIC(K)=QII(TOU(K)+(X(K)-TH*RK)*PK*(1-B*LOG(PK/PC)*M))
!    ---------------------------------------------------------------
!    IN  K      :  PLAN DE PROJECTION K = 1, 2 OU 3
!        MATER  :  COEFFICIENTS MATERIAU
!        SIG    :  CHAMPS DE CONTRAINTES A T
!        VIN    :  VARIABLES INTERNES A T
!    OUT SEUIL  :  SEUIL DU MECANISME DEVIATOIRE
!   ------------------------------------------------------------------
#include "asterf_types.h"
#include "asterfort/infniv.h"
    integer(kind=8) :: k, ndt, ndi
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: i, j
    real(kind=8) :: mater(22, 2), sig(6), vin(*), seuil
    real(kind=8) :: un, r, epsvp, pcr, pa, tole
    real(kind=8) :: degr, beta, b, m, phi, pcref, ptrac
    real(kind=8) :: p, q, x(2)
    real(kind=8) :: tou(3), th(2), touc(2)
    real(kind=8) :: d12, dd, deux
    aster_logical :: debug
    parameter(un=1.d0)
    parameter(tole=1.d-7)
    parameter(degr=0.0174532925199d0)
!       ------------------------------------------------------------
    common/tdim/ndt, ndi
    common/meshuj/debug
!
    data d12, deux/0.5d0, 2.d0/
!
    call infniv(ifm, niv)
!
!
! ==================================================================
! --- VARIABLES INTERNES -------------------------------------------
! ==================================================================
    epsvp = vin(23)
    r = vin(4+k)
    x(1) = vin(5+4*k)
    x(2) = vin(6+4*k)
    th(1) = vin(7+4*k)
    th(2) = vin(8+4*k)
!
!
! ==================================================================
! --- CARACTERISTIQUES MATERIAU ------------------------------------
! ==================================================================
    beta = mater(2, 2)
    b = mater(4, 2)
    phi = mater(5, 2)
    pcref = mater(7, 2)
    pa = mater(8, 2)
    pcr = pcref*exp(-beta*epsvp)
    ptrac = mater(21, 2)
    m = sin(degr*phi)
!
!
! ==================================================================
! --- PROJECTION DANS LE PLAN DEVIATEUR K --------------------------
! ==================================================================
    j = 1
    do i = 1, ndi
        if (i .ne. k) then
            tou(j) = sig(i)
            j = j+1
        end if
    end do
!
    tou(3) = sig(ndt+1-k)
    dd = d12*(tou(1)-tou(2))
!
!
! ==================================================================
! --- CALCUL DE PK, QCK --------------------------------------------
! ==================================================================
    p = d12*(tou(1)+tou(2))
    p = p-ptrac
    tou(1) = dd
    tou(2) = -dd
!
    if ((p/pa) .le. tole) then
        if (debug) write (ifm, '(A)') 'HUJCDC :: LOG(P/PA) NON DEFINI'
        seuil = -1.d0
        goto 999
    end if
!
    touc(2) = tou(3)-(x(2)-r*th(2))*p*(un-b*log(p/pcr))*m
    touc(1) = tou(1)-(x(1)-r*th(1))*p*(un-b*log(p/pcr))*m
!
    q = touc(1)**deux+(touc(2)**deux)/deux
    q = sqrt(q)
!
! ==================================================================
! --- CALCUL DU SEUIL DU MECANISME CYCLIQUE DEVIATOIRE K -----------
! ==================================================================
    seuil = -q/m/p-r*(un-b*log(p/pcr))
!
999 continue
end subroutine
