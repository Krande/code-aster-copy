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
subroutine dpmata(mod, mater, alpha, dp, dpdeno, &
                  pplus, se, seq, plas, dsde)
    implicit none
#include "asterc/r8prem.h"
#include "asterfort/lcopli.h"
#include "asterfort/lcprte.h"
#include "asterfort/utmess.h"
    real(kind=8) :: mater(5, 2), dp, dpdeno, se(6), seq, dsde(6, 6)
    real(kind=8) :: plas, alpha, pplus
    character(len=8) :: mod
! =====================================================================
! --- MISE A JOUR DES CONTRAINTES -------------------------------------
! =====================================================================
    integer(kind=8) :: ii, jj, ndt, ndi
    real(kind=8) :: un, deux, trois, young, nu, troisk, deuxmu, dsede(6, 6)
    real(kind=8) :: bidon(6, 6), pmat1(6, 6), pmat2(6, 6), pmat3(6, 6), param1
    real(kind=8) :: pmat4(6, 6), vunite(6), vect1(6), vect2(6), vect3(6)
    real(kind=8) :: pult, quatre, neuf, mater2(5, 2)
    parameter(neuf=9.0d0)
    parameter(quatre=4.0d0)
    parameter(trois=3.0d0)
    parameter(deux=2.0d0)
    parameter(un=1.0d0)
! =====================================================================
    common/tdim/ndt, ndi
! =====================================================================
! --- AFFECTATION DES VARIABLES ---------------------------------------
! =====================================================================
    young = mater(1, 1)
    nu = mater(2, 1)
    troisk = young/(un-deux*nu)
    deuxmu = young/(un+nu)
    pult = mater(4, 2)
    dsde(:, :) = 0.0d0
! =====================================================================
! --- CAS ELASTIQUE ---------------------------------------------------
! =====================================================================
    if (plas .eq. 0.0d0) then
        call lcopli('ISOTROPE', mod, mater(1, 1), dsde)
        goto 999
    else
        if (plas .ne. 2.0d0 .or. pplus .lt. pult) then
! =====================================================================
! --- INITIALISATIONS DE MATRICES ET VECTEURS UTILES ------------------
! =====================================================================
            dsede(:, :) = 0.0d0
            bidon(:, :) = 0.0d0
            pmat1(:, :) = 0.0d0
            pmat2(:, :) = 0.0d0
            pmat3(:, :) = 0.0d0
            pmat4(:, :) = 0.0d0
            vunite(:) = 0.0d0
            vect1(:) = 0.0d0
            vect2(:) = 0.0d0
            vect3(:) = 0.0d0
! =====================================================================
! --- CALCUL DU VECTEUR UNITE -----------------------------------------
! =====================================================================
            do ii = 1, ndi
                vunite(ii) = un
            end do
            if (plas .eq. 1.0d0) then
! =====================================================================
! --- CAS PLASTIQUE ---------------------------------------------------
! =====================================================================
! --- CALCUL DE DSEDE -------------------------------------------------
! =====================================================================
                do ii = 1, ndi
                    do jj = 1, ndi
                        dsede(ii, jj) = -deuxmu/trois
                    end do
                end do
                do ii = 1, ndt
                    dsede(ii, ii) = dsede(ii, ii)+deuxmu
                end do
! =====================================================================
! --- CALCUL DE PMAT1 -------------------------------------------------
! =====================================================================
                if (seq .gt. r8prem()) then
                    param1 = un-trois*deuxmu*dp/deux/seq
                    pmat1(1:ndt, 1:ndt) = param1*dsede(1:ndt, 1:ndt)
                else
                    call utmess('F', 'ALGORITH3_38')
                end if
! =====================================================================
! --- CALCUL DE PMAT2 -------------------------------------------------
! =====================================================================
                param1 = troisk/trois
                call lcprte(vunite, vunite, bidon)
                pmat2(1:ndt, 1:ndt) = param1*bidon(1:ndt, 1:ndt)
! =====================================================================
! --- CALCUL DE PMAT3 -------------------------------------------------
! =====================================================================
                param1 = neuf*deuxmu*deuxmu*dp/quatre/seq/seq/seq
                call lcprte(se, se, bidon)
                pmat3(1:ndt, 1:ndt) = param1*bidon(1:ndt, 1:ndt)
! =====================================================================
! --- CALCUL DE PMAT4 -------------------------------------------------
! =====================================================================
                param1 = trois*deuxmu/deux/seq
                vect1(1:ndt) = param1*se(1:ndt)
                param1 = troisk*alpha
                vect2(1:ndt) = param1*vunite(1:ndt)
                vect3(1:ndt) = vect1(1:ndt)+vect2(1:ndt)
                param1 = -un/dpdeno
                call lcprte(vect3, vect3, bidon)
                pmat4(1:ndt, 1:ndt) = param1*bidon(1:ndt, 1:ndt)
! =====================================================================
! --- CALCUL DE L OPERATEUR TANGENT -----------------------------------
! =====================================================================
                bidon(1:ndt, 1:ndt) = pmat1(1:ndt, 1:ndt)+pmat2(1:ndt, 1:ndt)
                pmat1(1:ndt, 1:ndt) = bidon(1:ndt, 1:ndt)+pmat3(1:ndt, 1:ndt)
                dsde(1:ndt, 1:ndt) = pmat1(1:ndt, 1:ndt)+pmat4(1:ndt, 1:ndt)
            else if (plas .eq. 2.0d0) then
! =====================================================================
! --- CAS DE LA PROJECTION AU SOMMET ----------------------------------
! =====================================================================
                param1 = troisk/trois-troisk*troisk*alpha*alpha/dpdeno
                call lcprte(vunite, vunite, bidon)
                dsde(1:ndt, 1:ndt) = param1*bidon(1:ndt, 1:ndt)
            end if
        else
! =====================================================================
! --- CAS DE LA PROJECTION AU SOMMET AVEC P > P_ULT -------------------
! --- DANS CE CAS ON PROPOSE DE CONSIDERER L'OPERATEUR TANGENT A UN ---
! --- FACTEUR MULTIPLICATIF PRES, QUE L'ON PREND ARBITRAIREMENT EGAL --
! --- A YOUNG/10E6 ----------------------------------------------------
! =====================================================================
            mater2(1, 1) = mater(1, 1)/1.0d6
            mater2(2, 1) = mater(2, 1)
            mater2(3, 1) = mater(3, 1)
            mater2(1, 2) = mater(1, 2)
            mater2(2, 2) = mater(2, 2)
            mater2(3, 2) = mater(3, 2)
            mater2(3, 2) = mater(3, 2)
            call lcopli('ISOTROPE', mod, mater2(1, 1), dsde)
        end if
    end if
! =====================================================================
999 continue
! =====================================================================
end subroutine
