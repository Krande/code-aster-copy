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
! aslint: disable=W0413
! comparaison aver r8gaem et -r8gaem uniquement

subroutine pielas(BEHinteg, ndim, npg, kpg, compor, &
                  typmod, mate, lgpg, vim, epsm, &
                  epsp, epsd, sigma, etamin, etamax, &
                  tau, copilo)
!
    use Behaviour_type
!
    use endo_loca_module, only: ELE_law => CONSTITUTIVE_LAW, ELE_Init => Init
    use endo_loca_module, only: ELE_PathFollowing => PathFollowing

    use vmis_isot_nl_module, only: VINL_law => CONSTITUTIVE_LAW, VINL_Init => Init
    use vmis_isot_nl_module, only: VINL_PathFollowing => PathFollowing
!
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8gaem.h"
#include "asterfort/pipedo.h"
#include "asterfort/pipedp.h"
#include "asterfort/pipeds.h"
#include "asterfort/pipepl.h"
#include "asterfort/utmess.h"
#include "blas/daxpy.h"
!
    type(Behaviour_Integ), intent(in) :: BEHinteg
    integer(kind=8) :: ndim, kpg, npg
    integer(kind=8) :: mate
    character(len=8) :: typmod(*)
    character(len=16) :: compor(*)
    integer(kind=8) :: lgpg
    real(kind=8) :: vim(lgpg, npg)
    real(kind=8) :: epsm(6), epsp(6), epsd(6)
    real(kind=8) :: copilo(5, npg)
    real(kind=8) :: etamin, etamax, tau, dka
    real(kind=8) :: sigma(6)
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (PILOTAGE - PRED_ELAS/DEFORMATION)
!
! PILOTAGE PAR PREDICTION ELASTIQUE
!
! ----------------------------------------------------------------------
!
!
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  NPG    : NOMBRE DE POINTS DE GAUSS
! IN  KPG    : NUMERO DU POINT DE GAUSS
! IN  TYPMOD : TYPE DE MODELISATION
! IN  MATE   : MATERIAU CODE
! IN  COMPOR : COMPORTEMENT
! IN  LGPG   : "LONGUEUR" DES VARIABLES INTERNES POUR 1 POINT DE GAUSS
!             CETTE LONGUEUR EST UN MAJORANT DU NBRE REEL DE VAR. INT.
! IN  VIM    : VARIABLES INTERNES EN T-
! IN  EPSM   : DEFORMATIONS AU TEMPS MOINS
! IN  EPSP   : CORRECTION DE DEFORMATIONS DUES AUX CHARGES FIXES
! IN  EPSD   : CORRECTION DE DEFORMATIONS DUES AUX CHARGES PILOTEES
! IN  SIGMA  : CONTRAINTES AVEC SQRT(2)
! IN  ETAMIN : BORNE INF. PILOTAGE
! IN  ETAMAX : BORNE SUP. PILOTAGE
! IN  TAU    : 2ND MEMBRE DE L'EQUATION F(ETA)=TAU
! OUT COPILO : COEFFICIENTS A0 ET A1 POUR CHAQUE POINT DE GAUSS
! ---------------------------------------------------------------------
    integer(kind=8) :: ndimsi, nsol, sgn(2)
    real(kind=8) :: sol(2), eps0(2*ndim), eps1(2*ndim)
    type(ELE_LAW) :: ELE_ldc
    type(VINL_LAW) :: VINL_ldc
    character(len=16) :: option
    blas_int :: b_incx, b_incy, b_n
! ---------------------------------------------------------------------
!
! --- INITIALISATIONS
!
    ndimsi = 2*ndim
    option = 'PILO_PRED_ELAS'
!
! --- CALCUL SUIVANT COMPORTEMENT
!
    if (compor(1) .eq. 'VMIS_ISOT_TRAC' .or. compor(1) .eq. 'VMIS_ISOT_LINE') then
        call pipepl(ndim, compor(1), typmod, tau, mate, &
                    sigma, vim(1, kpg), epsp, epsd, copilo(1, kpg), &
                    copilo(2, kpg), copilo(3, kpg), copilo(4, kpg), copilo(5, kpg))
!
    else if (compor(1) .eq. 'ENDO_ISOT_BETON') then
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, epsm, b_incx, epsp, &
                   b_incy)
!
        if (etamin .eq. -r8gaem() .or. etamax .eq. r8gaem()) then
            call utmess('F', 'MECANONLINE_60', sk=compor(1))
        end if
!
        call pipeds(ndim, typmod, tau, mate, vim(1, kpg), &
                    epsm, epsp, epsd, etamin, etamax, &
                    copilo(1, kpg), copilo(2, kpg), copilo(3, kpg), copilo(4, kpg), &
                    copilo(5, kpg))
!
!
    else if (compor(1) .eq. 'ENDO_LOCA_EXP' .or. compor(1) .eq. 'VMIS_ISOT_NL') then
!
        eps0 = epsm(1:ndimsi)+epsp(1:ndimsi)
        eps1 = epsd(1:ndimsi)
!
        if (compor(1) .eq. 'ENDO_LOCA_EXP') then
            ELE_ldc = ELE_Init(ndimsi, option, 'NONE', kpg, 1, mate, 100, 0.d0, 0.d0)
            call ELE_PathFollowing(ELE_ldc, vim(1, kpg)+tau, eps0, eps1, etamin, &
                                   etamax, 1.d-6, nsol, sol, sgn)
            if (ELE_ldc%exception .ne. 0) call utmess('F', 'PILOTAGE_83')

        else if (compor(1) .eq. 'VMIS_ISOT_NL') then
            VINL_ldc = VINL_Init(ndimsi, option, 'NONE', kpg, 1, mate, 100, 1.d-6)
            ! Normalisation tant que le meme tau s'applique a tout le maillage
            dka = tau*5.d-2
            call VINL_PathFollowing(VINL_ldc, dka, vim(:, kpg), eps0, eps1, nsol, sol, sgn)
            if (VINL_ldc%exception .ne. 0) call utmess('F', 'PILOTAGE_83')
        end if
!
        if (nsol .eq. 0) then
            copilo(5, kpg) = 0.d0
        else if (nsol .eq. 1) then
            copilo(1, kpg) = tau-sgn(1)*sol(1)
            copilo(2, kpg) = sgn(1)
        else if (nsol .eq. 2) then
            copilo(1, kpg) = tau-sgn(1)*sol(1)
            copilo(2, kpg) = sgn(1)
            copilo(3, kpg) = tau-sgn(2)*sol(2)
            copilo(4, kpg) = sgn(2)
        end if
!
!
    else if (compor(1) .eq. 'ENDO_ORTH_BETON') then
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, epsm, b_incx, epsp, &
                   b_incy)
!
        if (etamin .eq. -r8gaem() .or. etamax .eq. r8gaem()) then
            call utmess('F', 'MECANONLINE_60', sk=compor(1))
        end if
!
        call pipedo(ndim, typmod, tau, mate, vim(1, kpg), &
                    epsm, epsp, epsd, etamin, etamax, &
                    copilo(1, kpg), copilo(2, kpg), copilo(3, kpg), copilo(4, kpg), &
                    copilo(5, kpg))
    else if (compor(1) .eq. 'BETON_DOUBLE_DP') then
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, epsm, b_incx, epsp, &
                   b_incy)
!
        call pipedp(BEHinteg, kpg, 1, ndim, typmod, &
                    mate, epsm, sigma, vim(1, kpg), epsp, &
                    epsd, copilo(1, kpg), copilo(2, kpg))
!
    else
        call utmess('F', 'PILOTAGE_88', sk=compor(1))
    end if
end subroutine
