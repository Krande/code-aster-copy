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
subroutine drfnew(devg, devgii, traceg, dfds, dfdg, &
                  mu, k, dfdl)
!
    implicit none
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "blas/ddot.h"
    real(kind=8) :: devg(6), devgii, traceg, dfds(6), dfdg, mu, k, dfdl
! --- BUT : CALCUL DE DF/D(DELTA_LAMBDA) POUR NEWTON -------------------
! ======================================================================
! IN  : N      : NOMBRE TOTAL DE COMPOSANTES DU TENSEUR ----------------
! --- : ND     : NOMBRE DE COMPOSANTES DIAGONALES DU TENSEUR -----------
! --- : DEVG   : DEVIATEUR DE G ----------------------------------------
! --- : DEVGII : NORME DU DEVIATEUR DE G -------------------------------
! --- : TRACEG : PREMIER INVARIANT DE G --------------------------------
! --- : DFDS   : DERIVEE DE F PAR RAPPORT AUX CONTRAINTES --------------
! --- : DFDG   : DERIVEE DE F PAR RAPPORT A GAMP -----------------------
! --- : MU     : PARAMETRE MATERIAU ------------------------------------
! --- : K      : PARAMETRE MATERIAU ------------------------------------
! OUT : DFDL   : DF/DLAMBDA = - DF/DSIG.(2*MU*DEV(G) + K*TRACE(G)*I)
! ------------ :                + DF/DGAMP*RAC(2/3)*GII
! ======================================================================
    integer(kind=8) :: ii, ndt, ndi
    real(kind=8) :: vect1(6), scal1, mun, deux, trois
    blas_int :: b_incx, b_incy, b_n
! ======================================================================
! --- INITIALISATION DE PARAMETRES -------------------------------------
! ======================================================================
    parameter(mun=-1.0d0)
    parameter(deux=2.0d0)
    parameter(trois=3.0d0)
! ======================================================================
    common/tdim/ndt, ndi
! ======================================================================
    call jemarq()
! ======================================================================
! --- CALCUL INTERMEDIAIRE ---------------------------------------------
! ======================================================================
    do ii = 1, ndt
        vect1(ii) = deux*mu*devg(ii)
    end do
    do ii = 1, ndi
        vect1(ii) = vect1(ii)+k*traceg
    end do
    b_n = to_blas_int(ndt)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    scal1 = ddot(b_n, dfds, b_incx, vect1, b_incy)
! ======================================================================
! --- CALCUL FINAL -----------------------------------------------------
! ======================================================================
    dfdl = mun*scal1+dfdg*sqrt(deux/trois)*devgii
! ======================================================================
    call jedema()
! ======================================================================
end subroutine
