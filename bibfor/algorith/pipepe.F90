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
! person_in_charge: mickael.abbas at edf.fr
! aslint: disable=W1504
!
subroutine pipepe(BEHinteg, pilo, ndim, nno, npg, &
                  ipoids, ivf, idfde, geom, typmod, &
                  mate, compor, lgpg, deplm, sigm, &
                  vim, ddepl, depl0, depl1, copilo, &
                  iborne, ictau)
!
    use Behaviour_type
!
    implicit none
!
#include "jeveux.h"
#include "asterc/matfpe.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/pidefo.h"
#include "asterfort/pielas.h"
#include "asterfort/pipdef.h"
#include "asterfort/r8inir.h"
#include "blas/dcopy.h"
!
    type(Behaviour_Integ), intent(in) :: BEHinteg
    integer(kind=8) :: ndim, nno, npg
    integer(kind=8) :: mate, ipoids, ivf, idfde
    integer(kind=8) :: lgpg, iborne, ictau
    character(len=8) :: typmod(*)
    character(len=16) :: pilo, compor(*)
    real(kind=8) :: geom(ndim, *), deplm(*), ddepl(*)
    real(kind=8) :: sigm(2*ndim, npg), vim(lgpg, npg)
    real(kind=8) :: depl0(*), depl1(*)
    real(kind=8) :: copilo(5, npg)
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (PILOTAGE)
!
! CALCUL  DES COEFFICIENTS DE PILOTAGE POUR PRED_ELAS/DEFORMATION
!
! ----------------------------------------------------------------------
!
!
! IN  PILO   : MODE DE PILOTAGE: DEFORMATION, PRED_ELAS
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  NNO    : NOMBRE DE NOEUDS DE L'ELEMENT
! IN  NPG    : NOMBRE DE POINTS DE GAUSS
! IN  IPOIDS : POIDS DES POINTS DE GAUSS
! IN  IVF    : VALEUR DES FONCTIONS DE FORME
! IN  IDFDE  : DERIVEE DES FONCTIONS DE FORME ELEMENT DE REFERENCE
! IN  GEOM   : COORDONEES DES NOEUDS
! IN  TYPMOD : TYPE DE MODELISATION
! IN  MATE   : MATERIAU CODE
! IN  COMPOR : COMPORTEMENT
! IN  LGPG   : "LONGUEUR" DES VARIABLES INTERNES POUR 1 POINT DE GAUSS
!             CETTE LONGUEUR EST UN MAJORANT DU NBRE REEL DE VAR. INT.
! IN  DEPLM  : DEPLACEMENT EN T-
! IN  DDEPL  : INCREMENT DE DEPLACEMENT A L'ITERATION NEWTON COURANTE
! IN  SIGM   : CONTRAINTES DE CAUCHY EN T-
! IN  VIM    : VARIABLES INTERNES EN T-
! IN  DEPL0  : CORRECTION DE DEPLACEMENT POUR FORCES FIXES
! IN  DEPL1  : CORRECTION DE DEPLACEMENT POUR FORCES PILOTEES
! IN  IBORNE : ADRESSE JEVEUX POUR BORNES PILOTAGE
! IN  ICTAU  : ADRESSE JEVEUX POUR PARAMETRE PILOTAGE
! OUT COPILO : COEFFICIENTS A0 ET A1 POUR CHAQUE POINT DE GAUSS
!
!
!
!
!
    integer(kind=8) :: kpg, k, ndimsi
    real(kind=8) :: fm(3, 3), epsm(6), epsp(6), epsd(6)
    real(kind=8) :: rac2
    real(kind=8) :: etamin, etamax, tau, sigma(6)
    real(kind=8) :: dfdi(27, 3)
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
    call matfpe(-1)
!
! --- INITIALISATIONS
!
    rac2 = sqrt(2.d0)
    ndimsi = 2*ndim
!
    call r8inir(6, 0.d0, sigma, 1)
    call r8inir(npg*5, r8vide(), copilo, 1)
!
! --- TRAITEMENT DE CHAQUE POINT DE GAUSS
!
    do kpg = 1, npg
!
! --- CALCUL DES DEFORMATIONS
!
        call pipdef(ndim, nno, kpg, ipoids, ivf, &
                    idfde, geom, typmod, compor, deplm, &
                    ddepl, depl0, depl1, dfdi, fm, &
                    epsm, epsp, epsd)
!
! --- PILOTAGE PAR L'INCREMENT DE DEFORMATION
!
        if (pilo .eq. 'DEFORMATION') then
!
            call pidefo(ndim, npg, kpg, compor, fm, &
                        epsm, epsp, epsd, copilo)
!
!
! --- PILOTAGE PAR LA PREDICTION ELASTIQUE
!
        else if (pilo .eq. 'PRED_ELAS') then
            tau = zr(ictau)
            etamin = zr(iborne+1)
            etamax = zr(iborne)
!
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, sigm(1, kpg), b_incx, sigma, b_incy)
            do k = 4, ndimsi
                sigma(k) = sigma(k)*rac2
            end do
!
            call pielas(BEHinteg, ndim, npg, kpg, compor, &
                        typmod, mate, lgpg, vim, epsm, &
                        epsp, epsd, sigma, etamin, etamax, &
                        tau, copilo)
        else
            ASSERT(.false.)
        end if
!
    end do
!
    call matfpe(1)
!
end subroutine
