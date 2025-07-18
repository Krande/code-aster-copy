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

subroutine mbpk2c(sens, alpha, beta, h, covaini, jacini, jacdef, sigin, sigout)
!
    implicit none
#include "asterfort/assert.h"
#include "jeveux.h"
#include "asterfort/promat.h"
#include "asterfort/r8inir.h"
#include "asterc/r8prem.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: sens
    real(kind=8) :: covaini(3, 3), jacini, jacdef
    real(kind=8) :: alpha, beta, h
    real(kind=8) :: sigin(3), sigout(3)
! ----------------------------------------------------------------------
!    PASSAGE DES CONTRAINTES DE PIOLA-KIRCHHOFF II DANS LE REPÈRE NATUREL (LOCAL)
!    AUX CONTRAINTES INTEGREES DE CAUCHY DANS LE REPÈRE GLOBAL, POUR LES MEMBRANES
! ----------------------------------------------------------------------
! IN  SENS              0: DE PK2 A CAUCHY ; 1: DE CAUCHY A PK2
!     H                 EPAISSEUR DE LA MEMBRANE
!     COVAINI           COORD DES VECTEURS DE LA BASE COVARIANTES (RESP. ETAT INITIAL)
!     JACINI, JACDEF    JACOBIEN DE LA METRIQUE COVARIANTE        (RESP. ETAT INITIAL/DEFORME)
!     SIGIN            CONTRAINTES DE PIOLA KIRCHHOFF II DANS LE REPERE NATUREL
!
! OUT SIGOUT            CONTRAINTES DE CAUCHY DANS LE REPERE GOLBAL (S11,S22,S12)
! ----------------------------------------------------------------------
!
    integer(kind=8) :: i, k
    real(kind=8) :: sigintemp(3, 3), sigintemp1(3, 3), sigintemp2(3, 3)
    real(kind=8) :: vdirec(3), vortho(3)
    real(kind=8) :: matpt(3, 3), matp(3, 3), matptinv(3, 3), matpinv(3, 3)
    real(kind=8) :: projn, denomi, det

    call r8inir(3*3, 0.d0, sigintemp, 1)
    sigintemp(1, 1) = sigin(1)
    sigintemp(2, 2) = sigin(2)
    sigintemp(1, 2) = sigin(3)
    sigintemp(2, 1) = sigin(3)
!
! - CALCUL ET PROJECTION DU VECTEUR DIRECTION (X) SUR LA SURFACE
!
    vdirec(1) = cos(beta)*cos(alpha)
    vdirec(2) = cos(beta)*sin(alpha)
    vdirec(3) = -sin(beta)
!
    projn = 0.d0
    do i = 1, 3
        projn = projn+vdirec(i)*covaini(i, 3)
    end do
!
    if (abs(1.d0-projn*projn) .le. r8prem()) then
        call utmess('F', 'ELEMENTS_3')
    end if
!
    denomi = sqrt(1.d0-projn*projn)
    do i = 1, 3
        vdirec(i) = (vdirec(i)-projn*covaini(i, 3))/denomi
    end do

!
! - CALCUL DU VECTEUR TANGENT ORTHOGONAL AU VECTEUR DIRECTION PROJETE
!
    vortho(1) = covaini(2, 3)*vdirec(3)-covaini(3, 3)*vdirec(2)
    vortho(2) = covaini(3, 3)*vdirec(1)-covaini(1, 3)*vdirec(3)
    vortho(3) = covaini(1, 3)*vdirec(2)-covaini(2, 3)*vdirec(1)

! - ON CALCUL LA MATRICE DE PASSAGE DE LA BASE NATURELLE A BASE UTILISATEUR
    call r8inir(3*3, 0.d0, matp, 1)

    do i = 1, 2
        do k = 1, 3
            matp(1, i) = matp(1, i)+vdirec(k)*covaini(k, i)
            matp(2, i) = matp(2, i)+vortho(k)*covaini(k, i)
        end do
    end do
    matp(3, 3) = 1

    matpt = transpose(matp)
!
! - DIFFERENCIATION PK2->CAUCHY ET CAUCHY->PK2
!
    if (sens .eq. 0) then
! -- MATRICE DE PASSAGE DE BASE NATURELLE A BASE UTILISATEUR
!
        call promat(sigintemp, 3, 3, 3, matpt, 3, 3, 3, sigintemp1)
        call promat(matp, 3, 3, 3, sigintemp1, 3, 3, 3, sigintemp2)

        sigout(1) = h*sigintemp2(1, 1)*jacini/jacdef
        sigout(2) = h*sigintemp2(2, 2)*jacini/jacdef
        sigout(3) = h*sigintemp2(1, 2)*jacini/jacdef

    elseif (sens .eq. 1) then
! -- MATRICE DE PASSAGE DE BASE PROJETE A LA BASE NATURELLE
!
        det = (matp(1, 1)*matp(2, 2)-matp(1, 2)*matp(2, 1))

        call r8inir(3*3, 0.d0, matpinv, 1)
        matpinv(1, 1) = matp(2, 2)/det
        matpinv(2, 2) = matp(1, 1)/det
        matpinv(1, 2) = -matp(1, 2)/det
        matpinv(2, 1) = -matp(2, 1)/det
        matpinv(3, 3) = 1/det

        matptinv = transpose(matpinv)

        call promat(sigintemp, 3, 3, 3, matptinv, 3, 3, 3, sigintemp1)
        call promat(matpinv, 3, 3, 3, sigintemp1, 3, 3, 3, sigintemp2)

        sigout(1) = sigintemp2(1, 1)*jacdef/(jacini*h)
        sigout(2) = sigintemp2(2, 2)*jacdef/(jacini*h)
        sigout(3) = sigintemp2(1, 2)*jacdef/(jacini*h)

    else
        ASSERT(.false.)
    end if

end subroutine
