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

subroutine piesgv(neps, tau, mat, lccrma, vim, &
                  epsm, epsp, epsd, typmod, lcesga, &
                  etamin, etamax, lcesbo, copilo)
    implicit none
#include "asterf_types.h"
#include "asterfort/lcesma.h"
#include "asterfort/lcesvf.h"
#include "asterfort/piesfg.h"
#include "asterfort/utmess.h"
    interface
        subroutine lccrma(mat, fami, kpg, ksp, poum)
            integer(kind=8), intent(in) :: mat, kpg, ksp
            character(len=1), intent(in) :: poum
            character(len=*), intent(in) :: fami
        end subroutine lccrma
!
        subroutine lcesga(mode, eps, gameps, dgamde, itemax, &
                          precvg, iret)
            integer(kind=8), intent(in) :: mode, itemax
            real(kind=8), intent(in) :: eps(6), precvg
            integer(kind=8), intent(out) :: iret
            real(kind=8), intent(out) :: gameps, dgamde(6)
        end subroutine lcesga
!
        subroutine lcesbo(ep0, ep1, l0, l1, etamin, &
                          etamax, vide, etam, etap)
            real(kind=8), intent(in) :: ep0(6), ep1(6), l0, l1, etamin, etamax
            aster_logical, intent(out) :: vide
            real(kind=8), intent(out) :: etam, etap
        end subroutine lcesbo
    end interface
!
    character(len=8), intent(in) :: typmod(*)
    integer(kind=8), intent(in) :: neps, mat
    real(kind=8), intent(in) :: tau, epsm(neps), epsd(neps), epsp(neps), etamin, etamax, vim(3)
    real(kind=8), intent(out) :: copilo(2, *)
! --------------------------------------------------------------------------------------------------
!     PILOTAGE PRED_ELAS POUR ENDO_SCALAIRE (EN GRAD_VARI)
! --------------------------------------------------------------------------------------------------
! IN  NEPS    DIMENSION DES DEFORMATIONS (3*NDIM+2)
! IN  MAT     NATURE DU MATERIAU
! IN  LCCRMA  ROUTINE POUR LECTURE DES PARAMETRES DU CRITERE
! IN  VIM     VARIABLES INTERNES EN T-
! IN  EPSM    CHAMP DE DEFORMATION EN T-
! IN  EPSP    INCREMENT FIXE
! IN  EPSD    INCREMENT PILOTE
! IN  LCESGA  ROUTINE POUR CALCUL DU CRITERE EN DEFORMATION
! IN  ETAMIN  BORNE INF DU PILOTAGE
! IN  ETAMAX  BORNE SUP DU PILOTAGE
! OUT COPILO  COEFFICIENTS DE PILOTAGE :
!               F := COPILO(1,1)+COPILO(2,1)*ETA = TAU
!               F := COPILO(1,2)+COPILO(2,2)*ETA = TAU
! ----------------------------------------------------------------------
!  ITEMAX: NOMBRE MAX D'ITERATIONS POUR LA METHODE DE NEWTON
!  ERRA  : ERREUR TOLEREE SUR A DANS LA LDC (CRIT CVG)
!  RED   : REDUCTION DE L'ERREUR POUR EN FAIRE UN CRITERE DE PRECISION
    integer(kind=8), parameter :: itemax = 100
    real(kind=8), parameter :: red = 1.d-2, erra = 1.d-6
! ----------------------------------------------------------------------
    aster_logical :: cplan, croiss, gauche, droite, vide
    integer(kind=8) :: ndim, ndimsi, i, n
    real(kind=8) :: coplan
    real(kind=8) :: etam, etap, etal, precvg, l0, l1, etm, etp
    real(kind=8) :: gm, dgm, gp, dgp, gl, dgl
! ----------------------------------------------------------------------
    real(kind=8) :: lambda, deuxmu, troisk, gamma, rigmin, pc, pr, epsth
    common/lcee/lambda, deuxmu, troisk, gamma, rigmin, pc, pr, epsth
! ----------------------------------------------------------------------
    real(kind=8) :: pk, pm, pp, pq
    common/lces/pk, pm, pp, pq
! ----------------------------------------------------------------------
    real(kind=8) :: ep0(6), ep1(6), phi0, phi1, a, drda, precga
    common/pies/ep0, ep1, phi0, phi1, a, drda, precga
! ----------------------------------------------------------------------
!
!
!
! -- INITIALISATION
!
    cplan = typmod(1) .eq. 'C_PLAN  '
    ndim = (neps-2)/3
    ndimsi = 2*ndim
    if (ndim .eq. 2) then
        ep0(5:6) = 0
        ep1(5:6) = 0
    end if
!
! -- LECTURE DES CARACTERISTIQUES MATERIAU
!
    call lcesma(mat, 'NONE', 1, 1, '+', &
                lccrma)
!
!
! -- ENDOMMAGEMENT CIBLE
!
    a = vim(1)+tau
!
!
!  NON PILOTABLE CAR TROP PRES DE L'ENDOMMAGEMENT ULTIME
    if (a .ge. 0.99) goto 999
    drda = lcesvf(1, a)
!
!
! -- EXTRACTION DES DEFORMATIONS ET DES PARAMETRES NON LOCAUX
!
    do i = 1, ndimsi
        ep0(i) = epsm(i)+epsp(i)
        ep1(i) = epsd(i)
    end do
!
    if (cplan) then
        coplan = -lambda/(lambda+deuxmu)
        ep0(3) = coplan*(ep0(1)+ep0(2))
        ep1(3) = coplan*(ep1(1)+ep1(2))
    end if
!
    phi0 = epsm(ndimsi+2)+pr*epsm(ndimsi+1)+epsp(ndimsi+2)+pr*epsp(ndimsi+1)
    phi1 = epsd(ndimsi+2)+pr*epsd(ndimsi+1)
!
!
!
! -- AFFINER LES BORNES
!
    l0 = pm/pk/drda*(pk+pr*a-phi0)
    l1 = -pm/pk/drda*phi1
    call lcesbo(ep0, ep1, l0, l1, etamin, &
                etamax, vide, etm, etp)
!
!  PAS DE SOLUTION POUR LE PILOTAGE
    if (vide) then
        copilo(1, 3) = 0
        goto 999
    end if
!
!
!
! -- ESTIMATION INITIALE DU CRITERE DE CONVERGENCE
!
    precvg = red*pr*erra
    precga = min(red*precvg/abs(drda), pk/pm*1.d-3)
!
!
! -- RESOLUTION EN FONCTION DES BORNES
!
    call piesfg(lcesga, etm, gm, dgm)
    call piesfg(lcesga, etp, gp, dgp)
!
!
! 1. NE CONTRIBUE PAS AU PILOTAGE (TOUJOURS SOUS LE SEUIL)
!
    if (gm .le. 0 .and. gp .le. 0) then
        goto 999
    end if
!
!
! 2. BORNES SUPERIEURES AU SEUIL : DOUBLE NEWTON
!
    if (gm .ge. 0 .and. gp .ge. 0) then
!
!          INITIALISATION A GAUCHE ET A DROITE
        etam = etm
        etap = etp
        do n = 1, itemax
!
!              TEST DE CONVERGENCE
            gauche = abs(gm) .le. precvg
            droite = abs(gp) .le. precvg
            if (gauche .and. droite) goto 150
!
!              ABSENCE DE SOLUTION SI FONCTION AU-DESSUS DE ZERO
            if (dgm .ge. 0 .or. dgp .le. 0) then
                copilo(1, 3) = 0
                goto 999
            end if
!
!              METHODE DE NEWTON A GAUCHE ET A DROITE
            if (.not. gauche) etam = etam-gm/dgm
            if (.not. droite) etap = etap-gp/dgp

!              ABSENCE DE SOLUTION SI FONCTION AU-DESSUS DE ZERO
            if (etap .lt. etam) then
                copilo(1, 3) = 0
                goto 999
            end if
!
!              CALCUL DE LA FONCTION ET DERIVEE
            if (.not. gauche) call piesfg(lcesga, etam, gm, dgm)
            if (.not. droite) call piesfg(lcesga, etap, gp, dgp)
!
        end do
!
!          ECHEC DE LA RESOLUTION AVEC LE NOMBRE D'ITERATIONS REQUIS
        call utmess('F', 'PILOTAGE_83')
!
!
!          POST-TRAITEMENT DES SOLUTIONS
150     continue
        copilo(1, 1) = tau+etam
        copilo(2, 1) = -1
        copilo(1, 2) = tau-etap
        copilo(2, 2) = 1
        goto 999
    end if
!
!
!     3. BORNES DE PART ET D'AUTRE DU SEUIL --> NEWTON DEPUIS POSITIVE
    if (gm .ge. 0) then
        croiss = .false.
        etal = etm
        gl = gm
        dgl = dgm
    else
        croiss = .true.
        etal = etp
        gl = gp
        dgl = dgp
    end if
!
    do n = 1, itemax
!
!          TEST DE CONVERGENCE
        if (abs(gl) .le. precvg) goto 250
!
!          METHODE DE NEWTON A GAUCHE ET A DROITE
        etal = etal-gl/dgl
        call piesfg(lcesga, etal, gl, dgl)
    end do
!
    call utmess('F', 'PILOTAGE_83')
!
!      POST-TRAITEMENT DES SOLUTIONS
250 continue
    if (croiss) then
        copilo(1, 1) = tau-etal
        copilo(2, 1) = 1
    else
        copilo(1, 1) = tau+etal
        copilo(2, 1) = -1
    end if
!
999 continue
end subroutine
