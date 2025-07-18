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

subroutine glrc_calc_eps33(lambda, deuxmu, alfmc, gmt, gmc, &
                           tr2d, da1, da2, eps33, de33d1, &
                           de33d2, ksi2d, dksi1, dksi2, cof1, &
                           q2d, emp, cof2, dq2d)
!
!
! BUT : CALCUL DE LA COMPOSANTE DE LA DEFORMATION EZZ ET DE SES
!       DERIVEES POUR LE MODELE GLRC_DM
!
! IN:
!       LAMBDA  : PARAMETRE D ELASTICITE
!       DEUXMU  : PARAMETRE D ELASTICITE
!       TR2D    : TRACE 2D = EXX + EYY
!       D1      : ENDOMMAGEMENT DE LA PLAQUE D'UN COTE
!       D2      : ET DE L AUTRE
!       GMT     : PARAMETRE GAMMA POUR LA MEMBRANE EN TRACTION
!       GMC     : PARAMETRE GAMMA POUR LA MEMBRANE EN COMPRESSION
! OUT:
!       EPS33  : COMPOSANTE DE LA DEFORMATION EZZ
!       DE33D1 : DERIVEE DE EZZ PAR RAPPORT A D1
!       DE33D2 : DERIVEE DE EZZ PAR RAPPORT A D2
!       KSI2D : FONCTION CARACTERISTIQUE D ENDOMMAGEMENT,
!                KSI = KSI(TR2D,D1,D2)
! INTERNE:
!       DKSI1  : DERIVEE DE KSI PAR RAPPORT A D1
!       DKSI2  : DERIVEE DE KSI PAR RAPPORT A D2
! ----------------------------------------------------------------------
    implicit none
    integer(kind=8) :: k
    real(kind=8) :: tr2d, da1, da2, gmt, gmc, eps33, de33d1, de33d2
    real(kind=8) :: ksi2d, dksi1, dksi2, lambda, deuxmu, emp(2)
    real(kind=8) :: alfmc, mu
    real(kind=8) :: g1, g2, dg1, dg2, cof1(2), q2d(2)
    real(kind=8) :: gtr2(2), gi(2, 2), dgtr2(2), dgi(2, 2), cof2(2), dq2d(2)
!
    mu = 0.5d0*deuxmu
!
    if (tr2d .gt. 0.0d0) then
        g1 = (1.0d0+da1*gmt)/(1.0d0+da1)
        g2 = (1.0d0+da2*gmt)/(1.0d0+da2)
        dg1 = (1.0d0-gmt)/(1.0d0+da1)**2
        dg2 = (1.0d0-gmt)/(1.0d0+da2)**2
    else
        g1 = (alfmc+da1*gmc)/(alfmc+da1)
        g2 = (alfmc+da2*gmc)/(alfmc+da2)
        dg1 = alfmc*(1.0d0-gmc)/(alfmc+da1)**2
        dg2 = alfmc*(1.0d0-gmc)/(alfmc+da2)**2
    end if
!
! --  CALCUL DE KSIM (=KSI2D) ET DE SES DERIVEES PAR RAPPORT A D1 ET D2
!
    ksi2d = (g1+g2)*0.5d0
!
    dksi1 = -0.5d0*dg1
    dksi2 = -0.5d0*dg2
!
! --  CALCUL DE EPS33 ET DE SES DERIVEES PAR RAPPORT A D1 ET D2
!
    eps33 = -(lambda*tr2d*ksi2d)/(deuxmu+lambda*ksi2d)
!
    de33d1 = -(lambda*(tr2d+eps33)*dksi1)/(deuxmu+lambda*ksi2d)
    de33d2 = -(lambda*(tr2d+eps33)*dksi2)/(deuxmu+lambda*ksi2d)
!
! -------- CALCUL DE COF1 ET Q2D -----------
!
    if (tr2d .gt. 0.0d0) then
        gtr2(1) = 1.0d0-gmt
        gtr2(2) = 1.0d0-gmt
        dgtr2(1) = 0.0d0
        dgtr2(2) = 0.0d0
    else
        gtr2(1) = alfmc*(1.0d0-gmc)*(1.0d0+da1)**2/(alfmc+da1)**2
        gtr2(2) = alfmc*(1.0d0-gmc)*(1.0d0+da2)**2/(alfmc+da2)**2
        dgtr2(1) = 2.d0*alfmc*(1.0d0-gmc)*(1.0d0+da1)*(alfmc-1.0d0)
        dgtr2(1) = dgtr2(1)/(alfmc+da1)**3
        dgtr2(2) = 2.d0*alfmc*(1.0d0-gmc)*(1.0d0+da2)*(alfmc-1.0d0)
        dgtr2(2) = dgtr2(2)/(alfmc+da2)**3
    end if
!
    do k = 1, 2
        if (emp(k) .gt. 0.0d0) then
            gi(1, k) = 1.0d0-gmt
            gi(2, k) = 1.0d0-gmt
            dgi(1, k) = 0.0d0
            dgi(2, k) = 0.0d0
        else
            gi(1, k) = alfmc*(1.0d0-gmc)*(1.0d0+da1)**2/(alfmc+da1)**2
            gi(2, k) = alfmc*(1.0d0-gmc)*(1.0d0+da2)**2/(alfmc+da2)**2
            dgi(1, k) = 2.d0*alfmc*(1.0d0-gmc)*(1.0d0+da1)*(alfmc-1.0d0)
            dgi(1, k) = dgi(1, k)/(alfmc+da1)**3
            dgi(2, k) = 2.d0*alfmc*(1.0d0-gmc)*(1.0d0+da2)*(alfmc-1.0d0)
            dgi(2, k) = dgi(2, k)/(alfmc+da2)**3
        end if
    end do
!
    cof1(1) = 0.5d0*lambda*gtr2(1)
    cof1(2) = 0.5d0*lambda*gtr2(2)
    cof2(1) = 0.5d0*lambda*dgtr2(1)
    cof2(2) = 0.5d0*lambda*dgtr2(2)
    q2d(1) = 0.5d0*mu*(emp(1)**2*gi(1, 1)+emp(2)**2*gi(1, 2))
    q2d(2) = 0.5d0*mu*(emp(1)**2*gi(2, 1)+emp(2)**2*gi(2, 2))
    dq2d(1) = 0.5d0*mu*(emp(1)**2*dgi(1, 1)+emp(2)**2*dgi(1, 2))
    dq2d(2) = 0.5d0*mu*(emp(1)**2*dgi(2, 1)+emp(2)**2*dgi(2, 2))
!
end subroutine
