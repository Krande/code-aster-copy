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
subroutine glrcdd(zimat, maxmp, minmp, matr, ep, &
                  surfgp, q, epst, deps, dsig, &
                  ecr, delas, dsidep, normm, normn, &
                  crit, codret)
    implicit none
!
!     ROUTINE DE POST ET PRE TRAITEMENT POUR L INTEGRATION DE LA LOI DE
!     COMPORTEMENT GLRC_DAMAGE
!
! IN  ZIMAT : NUMERO DU MATERIAU
! IN  MAXMP : MOMENTS LIMITES ELASTIQUES EN FLEXION POSITIVE
! IN  MINMP : MOMENTS LIMITES ELASTIQUES EN FLEXION NEGATIVE
! IN  MATR : TABLEAU DES PARAMETRES MATERIAUX HOMOGENEISES
! IN  EP : EPAISSEUR TOTALE
! IN  SURFGP : SURFACE ASSOCIEE AU POINT DE GAUSS
! IN  Q : MATRICE DE PASSAGE ORTHO -> LOCAL
! IN  EPST : TENSEUR DEFORMATION TOTALE
! IN  DEPS : INCREMENT DE DEFORMATION
! IN  NORMM : NORME SUR LA FONCTION MP = F(N)
! IN  NORMN : NORME SUR LA FONCTION MP = F(N)
!
! OUT DSIG : INCREMENT DE CONTRAINTE
! OUT ECR : TABLEAU DES VARIABLES INTERNES
! OUT DELAS : MATRICE ELASTIQUE
! OUT DSIDEP : MATRICE TANGENTE
! person_in_charge: sebastien.fayolle at edf.fr
!
#include "asterfort/glrcad.h"
#include "asterfort/multsy.h"
#include "asterfort/ortloc.h"
#include "asterfort/r8inir.h"
    integer(kind=8) :: i, j, zimat, codret
!
    real(kind=8) :: matr(*), maxmp(*), minmp(*), normm, normn
    real(kind=8) :: ep, surfgp, q(2, 2), deps(*), epst(*)
    real(kind=8) :: dsig(*), ecr(*), dsidep(6, *), delas(6, *)
    real(kind=8) :: t1ve(9), c1(6, 6), c2(6, 6)
    real(kind=8) :: tq(2, 2), depslo(6), depsor(6), dpspor(6)
    real(kind=8) :: dfloc(6), dforth(6), nbacor(6), ep2s6
    real(kind=8) :: curloc(3), curort(3), ddiss
    real(kind=8) :: alpha, beta, gamma, dmax1, dmax2, k1, k2, e, nu, qp1, qp2
    real(kind=8) :: mf1, mf2, dam1, dam2, cuvcup(3)
    real(kind=8) :: rpara(5), crit(*)
!
!      SAVE  C1, C2, ALPHA, BETA, GAMMA, DMAX1, DMAX2, K1, K2
!
    e = matr(6)
    nu = matr(7)
    mf1 = matr(8)
    mf2 = matr(9)
    qp1 = matr(10)
    qp2 = matr(11)
!
    alpha = nu*e*ep**3/24/(1-nu**2)
    beta = e*ep**3/24.0d0/(1.0d0+nu)
    gamma = matr(12)
!
    if (mf1 .gt. 0.d0) then
        dmax1 = (1.d0-qp1)/(qp1-gamma)
        dmax2 = (1.d0-qp2)/(qp2-gamma)
        k1 = (1.d0-gamma)/4*mf1**2/(alpha+beta)
        k2 = (1.d0-gamma)/4*mf2**2/(alpha+beta)
    else
        dmax1 = (1.d0-qp2)/(qp2-gamma)
        dmax2 = (1.d0-qp1)/(qp1-gamma)
        k1 = (1.d0-gamma)/4*mf2**2/(alpha+beta)
        k2 = (1.d0-gamma)/4*mf1**2/(alpha+beta)
    end if
!
    call r8inir(6*6, 0.0d0, c1, 1)
    call r8inir(6*6, 0.0d0, c2, 1)
!
    do i = 1, 6
        c1(i, i) = matr(15+i)
        c2(i, i) = matr(21+i)
    end do
!
    do j = 1, 2
        do i = 1, 2
            tq(i, j) = q(j, i)
        end do
    end do
!
!     PASSSAGE DES TAUX DE DEFORMATION DANS LE REPERE ORTHOTROPE
!
    do j = 1, 6
        depslo(j) = deps(j)
    end do
!
!   ATTENTION: DEPS(3)    = 2*EPS_XY
!   ATTENTION: DEPS(6)    = 2*KAPPA_XY
!
!     CALCUL DU TENSEUR DES TAUX DE DEFORMATION DANS LE REPERE D ORTHO
    depslo(3) = deps(3)/2.d0
    depslo(6) = deps(6)/2.d0
!
    call multsy(tq, depslo, q, depsor)
    call multsy(tq, depslo(4), q, depsor(4))
!
    depsor(3) = 2.d0*depsor(3)
    depsor(6) = 2.d0*depsor(6)
!
!     CALCUL DU TENSEUR DE COURBURE DANS LE REPERE D ORTHO
    do i = 1, 3
        curloc(i) = epst(i+3)
    end do
!
    curloc(3) = epst(6)/2.d0
!
    call multsy(tq, curloc, q, curort)
    curort(3) = 2.d0*curort(3)
!
!     bending : elastoplasticity with kinematic softening
!
    ep2s6 = ep*ep/6.d0
!
    do j = 1, 6
        nbacor(j) = ecr(j+13)
    end do
!
!       NBACOR(1:3) = ep * multsym(tq,sig(1),q)
!       NBACOR(4:6) = ecr(17:19)    ! internal variable
!
!      NBACOR  = | n - backn |
!                     | m - backm |
!
    dam1 = ecr(8)*dmax1
    dam2 = ecr(9)*dmax2
!
    do j = 1, 3
        cuvcup(j) = curort(j)-ecr(j+3)
    end do
!
    rpara(1) = alpha
    rpara(2) = beta
    rpara(3) = gamma
    rpara(4) = k1
    rpara(5) = k2
!
!     INTEGRATION DE LA LOI DE COMPORTEMENT
!     POUR L INCREMENT DE DEF DANS LE REPERE D ORTHO DEPSOR
    call glrcad(zimat, maxmp, minmp, delas, rpara, &
                dmax1, dmax2, dam1, dam2, cuvcup, &
                c1, c2, nbacor, depsor, dpspor, &
                dforth, ddiss, dsidep, normm, normn, &
                crit, codret)
!
!     CONSTRUCTION DE LA MATRICE DE PASSAGE POUR LA MATRICE TANGENTE
    t1ve(1) = q(1, 1)*q(1, 1)
    t1ve(4) = q(2, 1)*q(2, 1)
    t1ve(7) = q(1, 1)*q(2, 1)
    t1ve(2) = t1ve(4)
    t1ve(5) = t1ve(1)
    t1ve(8) = -t1ve(7)
    t1ve(3) = -t1ve(7)-t1ve(7)
    t1ve(6) = t1ve(7)+t1ve(7)
    t1ve(9) = t1ve(1)-t1ve(4)
!
!     PASSAGE DE LA MATRICE TANGENTE DU REPERE D'ORTHOTROPIE AU REPERE
!     LOCAL DE L'ELEMENT. CE PASSAGE SE FAIT PAR SOUS MATRICES
    call ortloc(dsidep, 0, 0, t1ve)
    call ortloc(dsidep, 3, 3, t1ve)
    call ortloc(dsidep, 0, 3, t1ve)
!
    do i = 1, 3
        do j = 4, 6
            dsidep(j, i) = dsidep(i, j)
        end do
    end do
!
!     TRANSFERT DANS LE REPERE LOCAL DES INCREMENT D EFFORT ORTHO
    call multsy(q, dforth, tq, dfloc)
    call multsy(q, dforth(4), tq, dfloc(4))
!
!     CALCUL DES INCREMENTS DE CONTRAINTE DANS L EPAISSEUR
    do j = 1, 3
        dsig(j) = dfloc(j)/ep
        dsig(j+3) = dfloc(j+3)/ep2s6
    end do
!
!     STOCKAGE DES DEFORMATIONS DANS LE REPERE ORTHOTROPE
!     CALCUL DU TENSEUR DES TAUX DE DEFORMATION DANS LE REPERE LOCAL
!      DPSPOR(3) = DPSPOR(3) / 2.D0
!      DPSPOR(6) = DPSPOR(6) / 2.D0
!
!      CALL MULTSY(Q,DPSPOR,TQ,DPSPLO)
!      CALL MULTSY(Q,DPSPOR(4),TQ,DPSPLO(4))
!
!      DPSPLO(3)= 2.D0*DPSPLO(3)
!      DPSPLO(6) = 2.D0* DPSPLO(6)
!
    do j = 1, 6
        ecr(j) = ecr(j)+dpspor(j)
!        ECR(J) = ECR(J) + DPSPLO(J)
    end do
!
!     STOCKAGE DE L ENERGIE DISSIPEE PAR LA PLASTICITE
    ecr(7) = ecr(7)+surfgp*ddiss
!
!     STOCKAGE DE L ENERGIE DISSIPEE PAR ENDOMMAGEMENT
    ecr(10) = ecr(10)+surfgp*(k1*(dam1-ecr(8)*dmax1)+k2*(dam2-ecr(9)*dmax2))
!
!     STOCKAGE DE L ENDOMMAGEMENT
    ecr(8) = dam1/dmax1
    ecr(9) = dam2/dmax2
!
!     STOCKAGE DES EFFORTS ET MOMENTS DE RAPPEL CINEMATIQUE
!      NBACOR(3) = NBACOR(3) / 2.D0
!      NBACOR(6) = NBACOR(6) / 2.D0
!
!      CALL MULTSY(Q,NBACOR,TQ,NBACLO)
!      CALL MULTSY(Q,NBACOR(4),TQ,NBACLO(4))
!
!      NBACLO(3)= 2.D0*NBACLO(3)
!      NBACLO(6) = 2.D0* NBACLO(6)
!
    do j = 1, 6
        ecr(j+13) = nbacor(j)
!        ECR(J+13) = NBACLO(J)
    end do
!
end subroutine
