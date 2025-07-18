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

subroutine hujpxd(k, mater, sig, vin, prox, &
                  proxc)
    implicit none
!    HUJEUX:  CRITERE DE PROXIMITE POUR LES SEUILS DEVIATOIRES
!                        MONOTONES ET CYCLIQUES
!             FD(K) = QII(K) + M*PK*RK*( 1 - B*LOG(PK/PC) )
!    ---------------------------------------------------------------
!    IN  K      :  PLAN DE PROJECTION CONSIDERE (K = 1 A 3)
!        SIG    :  TENSEUR DES CONTRAINTES
!        MATER  :  PARAMETRES MATERIAU
!        VIN    :  VARIABLES INTERNES = ( R, X )
!    OUT PROX   :  CRITERE DE PROXIMITE
!                  .TRUE. MECANISMES ASSEZ PROCHES POUR ACTIVER
!                  LE MECANISME MONOTONE
!    ---------------------------------------------------------------
#include "asterf_types.h"
#include "asterfort/hujprj.h"
#include "asterfort/infniv.h"
    integer(kind=8) :: k, ndt, ndi
    integer(kind=8) :: ifm, niv
    real(kind=8) :: mater(22, 2), sig(6), vin(*)
    real(kind=8) :: un, r, epsvp, pcr, pa, tole1, tole2
    real(kind=8) :: degr, beta, b, m, phi, pcref, ptrac
    real(kind=8) :: sigd(3), p, q, dist, rh
    aster_logical :: debug, prox, proxc
    parameter(un=1.d0)
    parameter(tole1=1.d-6)
    parameter(tole2=1.d-7)
    parameter(degr=0.0174532925199d0)
!
    common/tdim/ndt, ndi
    common/meshuj/debug
!
    call infniv(ifm, niv)
!
!
! ==================================================================
! --- VARIABLES INTERNES -------------------------------------------
! ==================================================================
    epsvp = vin(23)
    rh = vin(k-4)
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
    ptrac = mater(21, 2)
    pcr = pcref*exp(-beta*epsvp)
    m = sin(degr*phi)
!
!
! ==================================================================
! --- PROJECTION DANS LE PLAN DEVIATEUR K ------------------------
! ==================================================================
    call hujprj(k-4, sig, sigd, p, q)
    if (((p-ptrac)/pa) .le. tole2) then
        if (debug) write (ifm, '(A)') 'HUJPXD :: LOG(P/PA) NON DEFINI'
        prox = .false.
        goto 999
    end if
!
!
! ==================================================================
! --- CALCUL DU SEUIL DU MECANISME DEVIATOIRE K ------------------
! ==================================================================
    r = -q/(m*(p-ptrac)*(un-b*log((p-ptrac)/pcr)))
    dist = abs(r-rh)/rh
!
    if (dist .lt. 1.d-5) then
        prox = .true.
    else
        prox = .false.
    end if
!
! ==================================================================
! --- SEUIL CYCLIQUE ELASTIQUE  + TANGENT AU SEUIL MONOTONE --------
! ==================================================================
    rh = sqrt(vin(4*k-11)**2+(vin(4*k-10)**2)/2.d0)
    if (rh .gt. tole1) then
        dist = abs(r-rh)/rh
    else
        dist = abs(r-rh)
    end if
!
    if ((dist .lt. tole1) .and. (vin(k) .eq. mater(18, 2))) then
        proxc = .true.
    end if
!
999 continue
end subroutine
