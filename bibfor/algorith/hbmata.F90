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
subroutine hbmata(se, dg, etap, i1e, sigeqe, &
                  vp, vecp, parame, derive, sig3, &
                  detadg, dgdl, nbmat, materf, dsidep)
    implicit none
#include "asterfort/calcdl.h"
#include "asterfort/lcprte.h"
    integer(kind=8) :: nbmat
    real(kind=8) :: se(6), dg, etap, i1e, dsidep(6, 6), materf(nbmat, 2)
    real(kind=8) :: vp(3), vecp(3, 3), sigeqe, parame(4), derive(5), sig3
    real(kind=8) :: detadg, dgdl
! ======================================================================
! -- HOEK BROWN : CALCUL DE LA MATRICE TANGENTE COHERENTE DSIG/DEPS ----
! ======================================================================
! IN  : NBMAT  : NOMBRE DE PARAMETRES MATERIAU -------------------------
! --- : MATERF : PARAMETRES MATERIAU -----------------------------------
! --- : SE     : DEVIATEUR DES CONTRAINTES ELASTIQUES ------------------
! --- : VP     : VALEURS PROPRES DU DEVIATEUR ELASTIQUE ----------------
! --- : VECP   : VECTEURS PROPRES DU DEVIATEUR ELASTIQUE ---------------
! --- : PARAME : VALEUR DES PARAMETRES DE LA LOI S*SIG, M*SIG, B -------
! --- : DERIVE : VALEUR DES DERIVEES DES PARAMETRES PAR RAPPORT A GAMMA
! --- : SIG3   : CONTRAINTE PRINCIPALE SIG3 ----------------------------
! --- : DG     : INCREMENT DU PARAMETRE D ECROUISSAGE GAMMA ------------
! --- : DETADG : DERIVEE DE ETA PAR RAPPORT A GAMMA --------------------
! --- : DGDL   : DERIVEE  DE GAMMA PAR RAPPORT A LAMBDA ----------------
! OUT : DSIDEP : DSIG/DEPS ---------------------------------------------
! ======================================================================
    integer(kind=8) :: ndt, ndi, ii, jj
    real(kind=8) :: un, deux, trois, mu, k
    real(kind=8) :: dsede(6, 6), param1, ddlde(6), seb(6)
    real(kind=8) :: vunite(6), bidon(6, 6), pmat1(6, 6), pmat6(6, 6)
    real(kind=8) :: pmat2(6, 6), pmat3(6, 6), pmat4(6, 6), pmat5(6, 6)
! ======================================================================
    parameter(deux=2.0d0)
    parameter(un=1.0d0)
    parameter(trois=3.0d0)
! ======================================================================
    common/tdim/ndt, ndi
! ======================================================================
! --- INITIALISATIONS --------------------------------------------------
! ======================================================================
    dsidep(:, :) = 0.0d0
    bidon(:, :) = 0.0d0
    dsede(:, :) = 0.0d0
    pmat1(:, :) = 0.0d0
    pmat2(:, :) = 0.0d0
    pmat3(:, :) = 0.0d0
    pmat4(:, :) = 0.0d0
    pmat5(:, :) = 0.0d0
    pmat6(:, :) = 0.0d0
    mu = materf(4, 1)
    k = materf(5, 1)
! =====================================================================
! --- CALCUL DU VECTEUR UNITE -----------------------------------------
! =====================================================================
    do ii = 1, ndi
        vunite(ii) = un
    end do
    do ii = ndi+1, 6
        vunite(ii) = 0.0d0
    end do
    do ii = 1, ndi
        seb(ii) = se(ii)
    end do
    do ii = ndi+1, ndt
        seb(ii) = se(ii)/sqrt(deux)
    end do
    do ii = ndt+1, 6
        seb(ii) = 0.0d0
    end do
! =====================================================================
! --- CALCUL DE DSEDE -------------------------------------------------
! =====================================================================
    do ii = 1, ndi
        do jj = 1, ndi
            dsede(ii, jj) = -deux*mu/trois
        end do
    end do
    do ii = 1, ndt
        dsede(ii, ii) = dsede(ii, ii)+deux*mu
    end do
! =====================================================================
! --- CALCUL DE K*DIEDE -----------------------------------------------
! =====================================================================
    call lcprte(vunite, vunite, bidon)
    pmat1(1:ndt, 1:ndt) = k*bidon(1:ndt, 1:ndt)
! =====================================================================
! --- CALCUL DE PARA*DSEDE --------------------------------------------
! =====================================================================
    param1 = un-trois*mu*dg/(sigeqe*(etap+un))
    pmat2(1:ndt, 1:ndt) = param1*dsede(1:ndt, 1:ndt)
    pmat6(1:ndt, 1:ndt) = pmat2(1:ndt, 1:ndt)+pmat1(1:ndt, 1:ndt)
! =====================================================================
! --- CALCUL DE SE*DSIGEQDE -------------------------------------------
! ====================================================================
    param1 = 9.0d0*mu*mu*dg/((etap+un)*sigeqe**3)
    call lcprte(seb, seb, bidon)
    pmat3 = param1*bidon
! ======================================================================
! --- CALCUL DE DDLAMBDA/DE ----------------------------------------
!=======================================================================
    call calcdl(vp, i1e, sigeqe, nbmat, materf, &
                parame, derive, sig3, vecp, etap, &
                dg, seb, detadg, dgdl, ddlde)
! ======================================================================
    param1 = trois*mu/sigeqe
    call lcprte(seb, ddlde, bidon)
    pmat4 = param1*bidon
! ======================================================================
    param1 = trois*k*(detadg*dgdl*dg/(etap+un)+etap)
    call lcprte(ddlde, vunite, bidon)
    pmat5 = param1*bidon
! ======================================================================
! --- CALCUL DE DSIG/DEPS ----------------------------------------------
! ======================================================================
    do ii = 1, ndt
        do jj = 1, ndt
            dsidep(ii, jj) = pmat6(ii, jj)+pmat3(ii, jj)-pmat4(ii, jj)-pmat5(ii, jj)
        end do
    end do
! ======================================================================
end subroutine
