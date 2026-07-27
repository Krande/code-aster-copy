! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
! aslint: disable=W1504
!
subroutine pipepe(BEHInteg, &
                  typmod, compor, &
                  pilo, ndim, nno, npg, &
                  ipoids, ivf, idfde, geom, &
                  lgpg, deplm, sigm, &
                  vim, ddepl, depl0, depl1, copilo, &
                  iborne, ictau)
!
    use Behaviour_type
    use Behaviour_module
    implicit none
!
#include "asterc/matfpe.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/pidefo.h"
#include "asterfort/pielas.h"
#include "asterfort/pipdef.h"
#include "asterfort/r8inir.h"
#include "blas/dcopy.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    character(len=8), intent(in) :: typmod(2)
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    integer(kind=8) :: ndim, nno, npg
    integer(kind=8) :: ipoids, ivf, idfde
    integer(kind=8) :: lgpg, iborne, ictau
    character(len=16) :: pilo
    real(kind=8) :: geom(ndim, *), deplm(*), ddepl(*)
    real(kind=8) :: sigm(2*ndim, npg), vim(lgpg, npg)
    real(kind=8) :: depl0(*), depl1(*)
    real(kind=8) :: copilo(5, npg)
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (PILOTAGE)
!
! CALCUL DES COEFFICIENTS DE PILOTAGE POUR PRED_ELAS/DEFORMATION
!
! --------------------------------------------------------------------------------------------------
!
! IN  PILO   : MODE DE PILOTAGE: DEFORMATION, PRED_ELAS
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  NNO    : NOMBRE DE NOEUDS DE L'ELEMENT
! IN  NPG    : NOMBRE DE POINTS DE GAUSS
! IN  IPOIDS : POIDS DES POINTS DE GAUSS
! IN  IVF    : VALEUR DES FONCTIONS DE FORME
! IN  IDFDE  : DERIVEE DES FONCTIONS DE FORME ELEMENT DE REFERENCE
! IN  GEOM   : COORDONEES DES NOEUDS
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
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    integer(kind=8), parameter :: ksp = 1
    integer(kind=8) :: kpg, k, ndimsi
    real(kind=8) :: fm(3, 3), epsm(6), epsp(6), epsd(6)
    real(kind=8) :: etamin, etamax, tau, sigma(6)
    real(kind=8) :: dfdi(MT_NNOMAX, 3)
    character(len=16) :: relaComp
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    call matfpe(-1)

! - INITIALISATIONS
    ndimsi = 2*ndim
    call r8inir(6, 0.d0, sigma, 1)
    call r8inir(npg*5, r8vide(), copilo, 1)
    relaComp = compor(RELA_NAME)

! - TRAITEMENT DE CHAQUE POINT DE GAUSS
    do kpg = 1, npg

! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(kpg, ksp, BEHInteg)

! ----- CALCUL DES DEFORMATIONS
        call pipdef(typmod, compor, &
                    ndim, nno, kpg, ipoids, ivf, &
                    idfde, geom, deplm, &
                    ddepl, depl0, depl1, dfdi, fm, &
                    epsm, epsp, epsd)

        if (pilo .eq. 'DEFORMATION') then
            call pidefo(compor, &
                        ndim, npg, kpg, fm, &
                        epsm, epsp, epsd, copilo)

        else if (pilo .eq. 'PRED_ELAS') then
            tau = zr(ictau)
            etamin = zr(iborne+1)
            etamax = zr(iborne)
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, sigm(1, kpg), b_incx, sigma, b_incy)
            do k = 4, ndimsi
                sigma(k) = sigma(k)*rac2
            end do
            call pielas(BEHInteg, &
                        typmod, relaComp, &
                        ndim, npg, kpg, &
                        lgpg, vim, epsm, &
                        epsp, epsd, sigma, etamin, etamax, &
                        tau, copilo)

        else
            ASSERT(ASTER_FALSE)
        end if
    end do
!
    call matfpe(1)
!
end subroutine
