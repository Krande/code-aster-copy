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
subroutine lcdpli(mod, nvi, option, materf, sigm, &
                  deps, vim, vip, sig, dsidep, &
                  iret)
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dplitg.h"
#include "asterfort/dpmata.h"
#include "asterfort/lcdevi.h"
#include "asterfort/lcopil.h"
#include "asterfort/lcopli.h"
#include "asterfort/majsig.h"
#include "asterfort/resdp1.h"
#include "asterfort/trace.h"
#include "blas/ddot.h"
    integer(kind=8) :: iret, nvi
    real(kind=8) :: deps(6), vim(nvi), vip(nvi), sig(6)
    real(kind=8) :: sigm(6), materf(5, 2), dsidep(6, 6)
    character(len=8) :: mod
    character(len=16) :: option
! =====================================================================
! --- LOI DE COMPORTEMENT DRUCKER PRAGER DE TYPE LINEAIRE -------------
! =====================================================================
! IN  OPTION  OPTION DE CALCUL (RAPH_MECA, RIGI_MECA_TANG OU FULL_MECA)
! IN  EPSM    CHAMP DE DEFORMATION EN T-
! IN  DEPS    INCREMENT DU CHAMP DE DEFORMATION
! IN  VIM     VARIABLES INTERNES EN T-
!               1   : ENDOMMAGEMENT (D)
!               2-7 : DEFORMATION ELASTIQUE (EPE AVEC RAC2 HD)
!               8   : INDICATEUR DISSIPATIF (1) OU ELASTIQUE (0)
! VAR VIP     VARIABLES INTERNES EN T+
!              IN  ESTIMATION (ITERATION PRECEDENTE)
!              OUT CALCULEES
! OUT SIGP    CONTRAINTES EN T+
! OUT DSIDEP  MATRICE TANGENTE
! OUT IRET    CODE RETOUR (0 = OK)
! =====================================================================
    aster_logical :: rigi, resi
    integer(kind=8) :: ndt, ndi, ii
    real(kind=8) :: trois, deux, dp, dpdeno, alpha, pmoins, pplus
    real(kind=8) :: hookf(6, 6), dkooh(6, 6), plas
    real(kind=8) :: epsp(6), epsm2(6), sige(6), se(6), siie, seq, i1e
    blas_int :: b_incx, b_incy, b_n
! =====================================================================
    parameter(deux=2.0d0)
    parameter(trois=3.0d0)
    common/tdim/ndt, ndi
! =====================================================================
! --- INITIALISATION --------------------------------------------------
! =====================================================================
    pmoins = vim(1)
    iret = 0
    resi = option(1:9) .eq. 'FULL_MECA' .or. option .eq. 'RAPH_MECA'
    rigi = option(1:9) .eq. 'FULL_MECA' .or. option(1:9) .eq. 'RIGI_MECA'
    ASSERT(resi .or. rigi)
! =====================================================================
! --- AFFECTATION DES VARIABLES ---------------------------------------
! =====================================================================
    alpha = materf(3, 2)
! =====================================================================
! --- OPERATEUR ELASTIQUE LINEAIRE ISOTROPE ---------------------------
! =====================================================================
    call lcopli('ISOTROPE', mod, materf(1, 1), hookf)
    call lcopil('ISOTROPE', mod, materf(1, 1), dkooh)
    epsm2(1:ndt) = matmul(dkooh(1:ndt, 1:ndt), sigm(1:ndt))
    epsp(1:ndt) = epsm2(1:ndt)+deps(1:ndt)
! =====================================================================
! --- INTEGRATION ELASTIQUE : SIGF = HOOKF EPSP -----------------------
! =====================================================================
    sige(1:ndt) = matmul(hookf(1:ndt, 1:ndt), epsp(1:ndt))
    call lcdevi(sige, se)
    b_n = to_blas_int(ndt)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    siie = ddot(b_n, se, b_incx, se, b_incy)
    seq = sqrt(trois*siie/deux)
    i1e = trace(ndi, sige)
!
! =====================================================================
! --- CALCUL DES CONTRAINTES ------------------------------------------
! =====================================================================
    if (resi) then
! =====================================================================
! --- RESOLUTION DU SYSTEME -------------------------------------------
! =====================================================================
        call resdp1(materf, seq, i1e, pmoins, dp, &
                    plas)
        if (plas .eq. 0.0d0) then
            do ii = 1, ndt
                sig(ii) = sige(ii)
            end do
        else
            call majsig(materf, se, seq, i1e, alpha, &
                        dp, plas, sig)
        end if
!
! =====================================================================
! --- STOCKAGE DES VARIABLES INTERNES ---------------------------------
! =====================================================================
        vip(1) = vim(1)+dp
        vip(2) = vim(2)+trois*alpha*dp
        vip(nvi) = plas
!
! =====================================================================
! --- PREPARATION AU CALCUL DE LA MATRICE TANGENTE --------------------
! =====================================================================
        dpdeno = dplitg(materf, vip(1), plas)
        pplus = vip(1)
    else
        plas = vim(nvi)
        dp = 0.0d0
        pplus = 0.0d0
        dpdeno = dplitg(materf, pmoins, plas)
    end if
! =====================================================================
! --- CALCUL DE LA MATRICE TANGENTE -----------------------------------
! =====================================================================
    if (rigi) then
        if (option(10:14) .eq. '_ELAS') then
            dsidep(1:ndt, 1:ndt) = hookf(1:ndt, 1:ndt)
        else
            call dpmata(mod, materf, alpha, dp, dpdeno, &
                        pplus, se, seq, plas, dsidep)
        end if
    end if
! =====================================================================
end subroutine
