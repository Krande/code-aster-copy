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

subroutine glrc_lc(epsm, deps, vim, option, sig, &
                   vip, dsidep, lambda, deuxmu, lamf, &
                   deumuf, gmt, gmc, gf, seuil, &
                   alf, alfmc, crit, &
                   epsic, epsiels, epsilim, codret, &
                   ep, is_param_opt, val_param_opt, t2iu)
! person_in_charge: sebastien.fayolle at edf.fr
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calc_glrcdm_err.h"
#include "asterfort/diago2.h"
#include "asterfort/glrc_calc_cst.h"
#include "asterfort/glrc_calc_eps33.h"
#include "asterfort/glrc_integ_loc.h"
#include "asterfort/glrc_sig_mat.h"
#include "asterfort/r8inir.h"
#include "asterfort/dxefro.h"
    integer(kind=8) :: codret
    real(kind=8) :: epsm(6), deps(6), vim(*), crit(*), seuil, alfmc
    real(kind=8) :: lambda, deuxmu, lamf, deumuf, alf, gmt, gmc, gf
    real(kind=8) :: epsic, epsiels, epsilim
    real(kind=8) :: sig(6), dsidep(6, 6), vip(*), vecp(2, 2), valp(2)
    character(len=16) :: option
    aster_logical :: is_param_opt(*), l_calc(2)
    real(kind=8) :: val_param_opt(*), ep, t2iu(4)
! ----------------------------------------------------------------------
!
!      LOI GLOBALE POUR LES PLAQUES/COQUES DKT - GLRC_DM
!
! IN:
!       LAMBDA  : PARAMETRE D ELASTICITE - MEMBRANE
!       DEUXMU  : PARAMETRE D ELASTICITE - MEMBRANE
!       LAMF    : PARAMETRE D ELASTICITE - FLEXION
!       DEUMUF  : PARAMETRE D ELASTICITE - FLEXION
!       GMT     : PARAMETRE GAMMA POUR LA MEMBRANE EN TRACTION
!       GMC     : PARAMETRE GAMMA POUR LA MEMBRANE EN COMPRESSION
!       GF      : PARAMETRE GAMMA POUR LA FLEXION
!       SEUIL   : INITIAL MEMBRANE
!       ALF     : PARAMETRE DE SEUIL FLEXION
!       EPSIC   : DEFORMATION AU PIC DE COMPRESSION
!       EPSIELS : DEFORMATION DES ACIERS A L'ETAT ULTIME DE SERVICE
!       EPSILIM : DEFORMATION A RUPTURE DES ACIERS
!       VIM     : VARIABLES INTERNES EN T-
!       OPTION  : TOUTES
!       CRIT    : CRITERES DE CONVERGENCE LOCAUX
!              (1) = NB ITERATIONS MAXI A CONVERGENCE
!                    (ITER_INTE_MAXI == ITECREL)
!              (2) = TYPE DE JACOBIEN A T+DT
!                    (TYPE_MATR_COMP == MACOMP)
!                     0 = EN VITESSE     >SYMETRIQUE
!                     1 = EN INCREMENTAL >NON-SYMETRIQUE
!              (3) = VALEUR TOLERANCE DE CONVERGENCE
!                    (RESI_INTE == RESCREL)
!              (5) = NOMBRE D'INCREMENTS POUR LE
!                    REDECOUPAGE LOCAL DU PAS DE TEMPS
!                    (ITER_INTE_PAS  == ITEDEC)
!                    -1,0,1 = PAS DE REDECOUPAGE
!                     N = NOMBRE DE PALIERS
!              (6) = TYPE D INTEGRATION LOCAL POUR LA LOI DE
!                    COMPORTEMENT (ALGO_INTE)
! OUT:
!       SIG     : CONTRAINTE
!       VIP     : VARIABLES INTERNES EN T+
!       DSIDEP  : MATRICE TANGENTE
!       CODRET  : CODE RETOUR DE L'INTEGRATION DE LA LDC
!                 0 => PAS DE PROBLEME
!                 1 => ABSENCE DE CONVERGENCE
! ----------------------------------------------------------------------
!
!       QM1 ET QM2 = Tm DANS R7.01.32
!       QFF        = Tf DANS R7.01.32
!
    aster_logical :: rigi, resi, coup
    aster_logical :: lelas, elas, elas1, elas2
    integer(kind=8) :: k, kdmax
    real(kind=8) :: eps(6), emp(2), efp(2), qff(2), eps8(8), epsu(6)
    real(kind=8) :: vmp(2, 2), vfp(2, 2), eps8out(8)
    real(kind=8) :: muf, trot, treps, eps33, de33d1, de33d2
    real(kind=8) :: da1, da2, ksi2d, dksi1, dksi2
    real(kind=8) :: tr2d, told, cof1(2), q2d(2)
    real(kind=8) :: cof2(2), dq2d(2), maxabs
    real(kind=8) :: tens(3), rx
!
! --  OPTION ET MODELISATION
    rigi = (option(1:4) .eq. 'RIGI' .or. option(1:4) .eq. 'FULL')
    resi = (option(1:4) .eq. 'RAPH' .or. option(1:4) .eq. 'FULL')
    coup = (option(6:9) .eq. 'COUP')
    if (coup) rigi = .true.
    lelas = option .eq. 'RIGI_MECA       '
!
! -- INITIALISATION
    if (lelas) then
        call r8inir(6, 0.d0, epsm, 1)
        call r8inir(6, 0.d0, deps, 1)
    end if
!
    muf = deumuf*0.5d0
!
! --  CALCUL DES EPSILON INITIAUX
    if (resi) then
        do k = 1, 6
            eps(k) = epsm(k)+deps(k)
        end do
    else
        do k = 1, 6
            eps(k) = epsm(k)
        end do
    end if
!
! --  ON UTILISE EPSILON SOUS FORME VECTORIELLE
! --  DONC ON DIVISE LES TERMES NON DIAGONNAUX PAR 2
    eps(3) = eps(3)/2.0d0
    eps(6) = eps(6)/2.0d0
!
! -- DIAGONALISATION DES DEFORMATIONS
    call diago2(eps(1), vmp, emp)
    call diago2(eps(4), vfp, efp)
!
! --  CALCUL DES TRACES
    tr2d = eps(1)+eps(2)
    trot = efp(1)+efp(2)
!
! --  CALCUL DES CONSTANTES INDEPENDANTES DE DA1, DA2 ET EPS33
    call glrc_calc_cst(lamf, muf, alf, gf, efp, qff)
!
! --  INITIALISATION DE DA1, DA2 ET EPS33
    if (lelas) then
        da1 = 0.0d0
        da2 = 0.0d0
    else
        da1 = vim(1)
        da2 = vim(2)
    end if
!
! --  EVOLUTION DE DA1, DA2 ET EPS33
!     INTEGRATION DE LA LOI DE COMPORTEMENT
    if (resi) then
        told = crit(3)
        kdmax = nint(crit(1))
!
        call glrc_integ_loc(lambda, deuxmu, seuil, alf, &
                            alfmc, gmt, gmc, cof1, &
                            vim, q2d, qff, tr2d, eps33, &
                            de33d1, de33d2, ksi2d, dksi1, dksi2, &
                            da1, da2, kdmax, told, codret, &
                            emp)
!
        if (da1 .lt. vim(1)) da1 = vim(1)
        if (da2 .lt. vim(2)) da2 = vim(2)
!
        elas1 = da1 .le. vim(1)
        elas2 = da2 .le. vim(2)
!
        elas1 = elas1 .or. lelas
        elas2 = elas2 .or. lelas
        elas = elas1 .and. elas2
!
        vip(1) = da1
        vip(2) = da2
        if (elas1) then
            vip(3) = 0.0d0
        else
            vip(3) = 1.0d0
        end if
        if (elas2) then
            vip(4) = 0.0d0
        else
            vip(4) = 1.0d0
        end if
        vip(5) = 1.d0-0.5d0*((1.d0+gmt*da1)/(1.d0+da1)+(1.d0+gmt*da2)/(1.d0+da2))
        vip(6) = 1.d0-0.5d0*((1.d0+gmc*da1)/(1.d0+da1)+(1.d0+gmc*da2)/(1.d0+da2))
        vip(7) = 1.d0-max((1.d0+gf*da1)/(1.d0+da1), (1.d0+gf*da2)/(1.d0+da2))

        if (is_param_opt(1)) then

!           passage des deformation dans le repere utilisateur
            eps8(1:6) = eps(1:6)
            eps8(7:8) = 0.d0
            call dxefro(1, t2iu, eps8, eps8out)
            epsu(1:6) = eps8out(1:6)

            rx = val_param_opt(1)

            maxabs = max(abs(epsu(1)-0.5d0*ep*rx*epsu(4)), &
                         abs(epsu(1)+0.5d0*ep*rx*epsu(4)))
            vip(8) = maxabs/epsiels
            vip(9) = maxabs/epsilim

            maxabs = max(abs(epsu(2)-0.5d0*ep*rx*epsu(5)), &
                         abs(epsu(2)+0.5d0*ep*rx*epsu(5)))
            vip(10) = maxabs/epsiels
            vip(11) = maxabs/epsilim

            tens(1) = epsu(1)-ep/2*epsu(4)
            tens(2) = epsu(2)-ep/2*epsu(5)
            tens(3) = epsu(3)-ep/2*epsu(6)
            call diago2(tens, vecp, valp)
            vip(12) = -min(valp(1), valp(2), 0.d0)/epsic

            tens(1) = epsu(1)+ep/2*epsu(4)
            tens(2) = epsu(2)+ep/2*epsu(5)
            tens(3) = epsu(3)+ep/2*epsu(6)
            call diago2(tens, vecp, valp)
            vip(13) = -min(valp(1), valp(2), 0.d0)/epsic

            vip(14) = max(vim(14), epsu(1), epsu(2), 0.d0)
            vip(15) = max(vim(15), -epsu(1), -epsu(2), 0.d0)
            vip(16) = max(vim(16), abs(epsu(4)), abs(epsu(5)))

            if (is_param_opt(2)) then

                if (vip(15) .gt. vim(15)) then
                    l_calc(1) = .true.
                else
                    l_calc(1) = .false.
                    vip(17) = vim(17)
                end if

                if (vip(16) .gt. vim(16)) then
                    l_calc(2) = .true.
                else
                    l_calc(2) = .false.
                    vip(18) = vim(18)
                end if

                call calc_glrcdm_err(l_calc, vip(15), vip(16), gf, &
                                     gmc, epsic, ep, val_param_opt, &
                                     vip(17), vip(18))

            else
                vip(17:18) = 0.d0
            end if
        else
            vip(8:18) = 0.d0
        end if
    else
        if (lelas) then
            da1 = 0.0d0
            da2 = 0.0d0
            elas1 = .true.
            elas2 = .true.
            elas = .true.
        else
            da1 = vim(1)
            da2 = vim(2)
            elas1 = nint(vim(3)) .eq. 0
            elas2 = nint(vim(4)) .eq. 0
            elas = (elas1 .and. elas2)
        end if
    end if
    call glrc_calc_eps33(lambda, deuxmu, alfmc, gmt, gmc, &
                         tr2d, da1, da2, eps33, de33d1, &
                         de33d2, ksi2d, dksi1, dksi2, cof1, &
                         q2d, emp, cof2, dq2d)
!
! --  CALCUL DE LA TRACE 3D
    treps = tr2d+eps33
!
! --  CALCUL DES CONTRAINTES GENERALISEES ET DE LA MATRICE TANGENTE
!
    call glrc_sig_mat(lambda, deuxmu, lamf, deumuf, alf, &
                      alfmc, emp, efp, eps, vmp, &
                      vfp, tr2d, trot, treps, gmt, &
                      gmc, gf, da1, da2, ksi2d, &
                      qff, cof1, q2d, de33d1, de33d2, &
                      elas, elas1, elas2, coup, rigi, &
                      resi, option, dsidep, sig, cof2, &
                      dq2d)
!
end subroutine
