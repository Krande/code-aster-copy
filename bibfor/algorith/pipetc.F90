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
subroutine pipetc(mat, sup, sud, mup, mud, &
                  vim, tau, copilo)
!
!
    implicit none
#include "asterc/r8vide.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "asterfort/zerop2.h"
    integer(kind=8) :: mat
    real(kind=8) :: sup(3), sud(3), mup(3), mud(3), tau, vim(*), copilo(5)
!
! ----------------------------------------------------------------------
!     PILOTAGE PRED_ELAS POUR LA LOI D'INTERFACE CZM_TAC_MIX
!
! IN  MAT    : MATERIAU
! IN  SUP    : SAUT DU AUX CHARGES FIXES
! IN  SUD    : SAUT DU AUX CHARGES PILOTEES
! IN  MUP    : MULTIPLICATEUR DU AUX CHARGES FIXES
! IN  MUD    : MULTIPLICATEUR DU AUX CHARGES PILOTEES
! IN  VIM    : VARIABLES INTERNES EN T-
! IN  TAU    : 2ND MEMBRE DE L'EQUATION F(ETA)=TAU
! OUT COPILO : COEFFICIENTS DU TIR ELASTIQUE LINEARISE AUTOUR DES SOL.
!                FEL = COPILO(1) + COPILO(2)*ETA
!                FEL = COPILO(3) + COPILO(4)*ETA
!                COPILO(5) <> R8VIDE => PAS DE SOLUTION
! ----------------------------------------------------------------------
    integer(kind=8) :: nrac, i
    real(kind=8) :: sc, gc, dc, h, r, ka, ga, sk, val(3), tmp
    real(kind=8) :: tt2, tp(3), td(3), tpn, tdn, c0, c1, c2, n0, n1, n2, rac(4)
    real(kind=8) :: eta, pente
    integer(kind=8) :: cod(3), kpg, spt
    character(len=16) :: nom(3)
    character(len=8) :: fami, poum
!
    data nom/'GC', 'SIGM_C', 'PENA_LAGR'/
! ----------------------------------------------------------------------
!
! -- SEUIL EN ENERGIE DISSIPE
    ga = vim(4)+tau
!
! -- CAS DE L'ENDOMMAGEMENT SATURE
    if (ga .gt. 1.d0) goto 999
!
! -- RECUPERATION DES PARAMETRES PHYSIQUES
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    call rcvalb(fami, kpg, spt, poum, mat, &
                ' ', 'RUPT_FRAG', 0, ' ', [0.d0], &
                3, nom, val, cod, 2)
    gc = val(1)
    sc = val(2)
    dc = 2.d0*gc/sc
    h = sc/dc
    r = h*val(3)
!
!    CALCUL DE KAPPA : KA = DC*(1-SQRT(1-GA)) : SEUIL EN SAUT
    tmp = sqrt(max(0.d0, 1.d0-ga))
    tmp = dc*(1.d0-tmp)
    tmp = max(0.d0, tmp)
    tmp = min(dc, tmp)
    ka = tmp
    sk = max(0.d0, sc-h*ka)
!
!   CALCUL DU SEUIL (EN CONTRAINTE)
    tt2 = (r*ka+sk)**2.d0
! -- CALCUL DU SECOND MEMBRE
!
!    FORCE COHESIVE AUGMENTEE
    tp(1) = mup(1)+r*sup(1)
    tp(2) = mup(2)+r*sup(2)
    tp(3) = mup(3)+r*sup(3)
!
!    FORCE COHESIVE AUGMENTEE
    td(1) = mud(1)+r*sud(1)
    td(2) = mud(2)+r*sud(2)
    td(3) = mud(3)+r*sud(3)
!
    tpn = tp(1)
    tdn = td(1)
!
!
! -- CALCUL DES FORMES QUADRATIQUES
!
!    CISAILLEMENT : C0 + 2.C1 ETA + C2 ETA**2
!    OUVERTURE    : N0 + 2.N1 ETA + N2 ETA**2
!
    c0 = tp(2)*tp(2)+tp(3)*tp(3)
    c1 = tp(2)*td(2)+tp(3)*td(3)
    c2 = td(2)*td(2)+td(3)*td(3)
!
    n0 = tpn*tpn
    n1 = tpn*tdn
    n2 = tdn*tdn
!
!    SI LE POINT N'EST PAS PILOTABLE :
    if (c2+n2 .eq. 0.d0) goto 999
!
!
! -- RESOLUTION DES EQUATIONS
!
!             C(ETA) = TT2
!    N(ETA) + C(ETA) = TT2
!
!    INITIALISATION
    call r8inir(4, r8vide(), rac, 1)
!
!    EQUATION SUR LE CISAILLEMENT SEUL
    if (c2 .ne. 0.d0) call zerop2(2.d0*c1/c2, (c0-tt2)/c2, rac(1), nrac)
!
!    EQUATION SUR LE CISAILLEMENT + L'OUVERTURE
    call zerop2(2.d0*(c1+n1)/(c2+n2), (c0+n0-tt2)/(c2+n2), rac(3), nrac)
!
!
! -- SELECTION DES SOLUTIONS CONFORMES AVEC LE SIGNE SUR L'OUVERTURE
!
    nrac = 0
    do i = 1, 4
        if (rac(i) .ne. r8vide()) then
            eta = rac(i)
!
            if (i .le. 2) then
                if ((tpn+eta*tdn) .lt. 0.d0) then
                    pente = 2.d0*(c2*eta+c1)
                    copilo(nrac+1) = tau-pente*eta
                    copilo(nrac+2) = pente
                    nrac = nrac+2
                end if
            else
                if ((tpn+eta*tdn) .ge. 0.d0) then
                    pente = 2.d0*((c2+n2)*eta+c1+n1)
                    copilo(nrac+1) = tau-pente*eta
                    copilo(nrac+2) = pente
                    nrac = nrac+2
                end if
            end if
!
        end if
    end do
    if (nrac .eq. 0) copilo(5) = 0.d0
!
999 continue
end subroutine
