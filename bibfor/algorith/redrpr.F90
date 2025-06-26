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
subroutine redrpr(mod, imate, sigp, vip, dsde, &
                  icode)
    implicit none
#include "asterfort/dplitg.h"
#include "asterfort/dpmata.h"
#include "asterfort/dpmate.h"
#include "asterfort/dppatg.h"
#include "asterfort/lcdevi.h"
#include "blas/ddot.h"
    integer(kind=8) :: imate, icode
    real(kind=8) :: vip(*), sigp(*), dsde(6, 6)
    character(len=8) :: mod
! =====================================================================
! --- ROUTINE POUR LES RECUPERATIONS DES DONNES -----------------------
! --- POUR LE CALCUL DU TENSEUR TANGENT -------------------------------
! --- ICODE = 0 CORRESPONDT AU CAS ELASTIQUE --------------------------
! --- ICODE = 1 SINON -------------------------------------------------
! =====================================================================
    integer(kind=8) :: ndt, ndi, nvi, typedp
    real(kind=8) :: pplus, materf(5, 2), hookf(6, 6), dpdeno, dp
    real(kind=8) :: se(6), seq, plas, alpha, phi
    real(kind=8) :: siie, deux, trois
    blas_int :: b_incx, b_incy, b_n
! ======================================================================
    common/tdim/ndt, ndi
! =====================================================================
! --- INITIALISATION --------------------------------------------------
! =====================================================================
    parameter(deux=2.0d0)
    parameter(trois=3.0d0)
! =====================================================================
    se(:) = 0.0d0
    dsde(:, :) = 0.0d0
    hookf(:, :) = 0.0d0
! =====================================================================
! --- RECUPERATION DES DONNEES MATERIAUX ------------------------------
! =====================================================================
    call dpmate(mod, imate, materf, ndt, ndi, &
                nvi, typedp)
    pplus = vip(1)
    dp = 0.0d0
    plas = vip(nvi)
    icode = 1
    seq = 0.0d0
    if (typedp .eq. 1) then
! =====================================================================
! --- RECUPERATION DE R' POUR UNE LOI DP DE TYPE LINEAIRE -------------
! =====================================================================
        alpha = materf(3, 2)
        dpdeno = dplitg(materf, pplus, plas)
    else if (typedp .eq. 2) then
! =====================================================================
! --- RECUPERATION DE R' POUR UNE LOI DP DE TYPE PARABOLIQUE ----------
! =====================================================================
        phi = materf(2, 2)
        alpha = deux*sin(phi)/(trois-sin(phi))
        dpdeno = dppatg(materf, pplus, plas)
    end if
    if (plas .eq. 0.0d0) then
        icode = 0
    else
        if (plas .eq. 1.0d0) then
! =====================================================================
! --- INTEGRATION ELASTIQUE : SIGF = HOOKF EPSP -----------------------
! =====================================================================
            call lcdevi(sigp, se)
            b_n = to_blas_int(ndt)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            siie = ddot(b_n, se, b_incx, se, b_incy)
            seq = sqrt(trois*siie/deux)
        end if
        call dpmata(mod, materf, alpha, dp, dpdeno, &
                    pplus, se, seq, plas, dsde)
    end if
! =====================================================================
! =====================================================================
end subroutine
